#!/usr/bin/env nextflow

/*
This workflow contains pre- and post-processing steps to 
calculate Magnetization Transfer Saturation Index (MTsat) map along
with a longitudinal relaxation time (T1) map.

Dependencies: 
    - Advanced notmarization tools (ANTs, https://github.com/ANTsX/ANTs)
        - Installation: Built from source 
    - qMRLab (https://qmrlab.org) 
        - MATLAB/Octave 
Docker: 
    - 


Author:
    Agah Karakuzu 2019
    agahkarakuzu@gmail.com 

Users: Please see USAGE for utilizatoin purposes.
 */

/*Set defaults for parameters determining logic flow to false*/
params.root = false 
params.bids = false 
params.help = false

/*Define bindings for --help*/
if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["test":"$params.dim",
                "test2":"$params.transform"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

/*Scrape file names from a BIDS-compatible dataset
Note:
    BIDS for qMRI is currently under development (BEP001,https://github.com/bids-standard/bep001)
    The current format is valid as of late 2019 and subjected to change.
    For B1plusmaps, there is not a specification yet. To circumvent this 
    issue, these (optional) maps are assumed to be located at the fmap
    folder with _B1plusmap suffix.   
*/
if(params.root && params.bids){
    log.info "Input: $params.root"
    root = file(params.root)
    /* Here, alphabetical indexes matter. Therefore, MToff -> MTon -> T1w */
    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_acq-{MToff,MTon,T1w}_MTS.nii.gz", maxDepth: 2, size: 3, flat: true)
    (pdw, mtw, t1w) = in_data
        .map{sid, MToff, MTon, T1w  -> [    tuple(sid, MToff),
                                            tuple(sid, MTon),
                                            tuple(sid, T1w)]}                                   
        .separate(3)
                 
    /* Look for B1map in fmap folders */
    Channel
           .fromPath("$root/**/fmap/sub-*_B1plusmap.nii.gz", maxDepth:2)   
           .map{[it.parent.parent.name, it]}
           .set{b1plus}      
             
}   

/*If data is not BIDS-compatible, define a custom directory structure
    Please see USAGE for further details. 
*/ 
else if(params.root && !params.bids){
    log.info "Input: $params.root"
    root = file(params.root)
    /* Here, alphabetical indexes matter. Therefore, MToff -> MTon -> T1w */
    in_data = Channel.fromFilePairs("$root/**/*{MTw.nii.gz,PDw.nii.gz,T1w.nii.gz}", maxDepth: 1, size: 3, flat: true){it.parent.name}

    (mtw, pdw, t1w) = in_data
        .map{sid, MTw, PDw, T1w  -> [    tuple(sid, MTw),
                                         tuple(sid, PDw),
                                         tuple(sid, T1w)]}                                   
        .separate(3)

    /* Look for B1map in fmap folders */
    Channel 
    .fromPath("$root/**/*B1plusmap.nii.gz",
                    maxDepth:1)                              
    .map{[it.parent.name, it]}
    .set{b1plus}


}
else{
    error "ERROR: Arguments (--root) and (--bids) must be passed. See USAGE."
}

/*Each data type is defined as a channel. To pass all the channels 
  to the same process accurately, these channels must be joined. 
*/ 

/*Split T1w into two channels*/
t1w.into{t1w_pre; t1w_post}

/* Merge PDw, MTw and T1w for alingment and brain extraction*/
pdw 
    .join(mtw)
    .join(t1w_pre)
    .set{mtsat_for_preproc}

log.info "qMRflow: MTsat pipeline"
log.info "======================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.info ""
log.info "DATA"
log.info "===="
log.info ""
if (params.bids){
log.info "BIDS option has been enabled."
log.warn "qMRI protocols will be read from sidecar .json files."
}
else{
log.info "Custom file/folder organization described in USAGE will be assumed."
log.warn "If an mtsat_protocol.json file is provided for a subject, acquisition metadata will be overridden. Otherwise, those defined in nextflow.config will be used."
}

log.info ""
log.info "OPTIONS"
log.info "======="
log.info ""
log.info "[ANTs Registration]"
log.info "-------------------"
log.info "Dimensionality: $params.dim"
log.info "Metric: $params.metric"
log.info "Weight: $params.metric_weight"
log.info "Number of bins: $params.metric_bins"
log.info "Sampling type: $params.metric_sampling"
log.info "Sampling percentage: $params.metric_samplingprct"
log.info "Transform: $params.transform"
log.info "Convergence: $params.convergence"
log.info "Shrink factors: $params.shrink"
log.info "Smoothing sigmas: $params.smoothing"
log.info ""
log.info "[qMRLab mt_sat]"
log.info "---------------"
log.info "NOTE: mt_sat protocol inputs are subjected to change based on the (--bids) option. "
log.info "\tSee more details at USAGE."
log.info ""
log.info "Default values provided in the nextflow.config: "
log.info ""
log.info "Flip angles:\n\t MTon: $params.mtw_fa\n\t MToff: $params.pdw_fa\n\t T1w: $params.t1w_fa"
log.info "Repetition times:\n\t MTon: $params.mtw_tr\n\t MToff: $params.pdw_tr\n\t T1w: $params.t1w_tr"
log.info ""
if (params.use_b1map){
log.info "B1map option has been enabled."  
log.warn "Process will be skipped for participants lacking a B1map."   
log.info "B1 correction factor: $params.b1_cor_factor"}
if (!params.use_b1map){
log.info "B1map option has been disabled."
log.warn "Process won't take any (possibly) existing B1maps into account."
}

/*Perform rigid registration to correct for head movement across scans:
    - MTw (moving) --> T1w (fixed)
    - PDw (moving) --> T1w (fixed)
*/     

process Align_And_Extract {
    cpus 2


    publishDir = "$root/derivatives/qMRLab"


    input:
    set sid, file(pdw), file(mtw), file(t1w) from mtsat_for_preproc

    output:
    set sid, "${sid}_acq-MTon-to-T1w_MTS.nii.gz", "${sid}_acq-MToff-to-T1w_MTS.nii.gz",\
    "${sid}_acq-T1w_mask.nii.gz" into mtsat_from_preproc

    script:
    """
    touch ${sid}_acq-MTon-to-T1w_MTS.nii.gz
    touch ${sid}_acq-MToff-to-T1w_MTS.nii.gz
    touch ${sid}_acq-T1w_mask.nii.gz
    """
}


/* Split t1w_post into two to deal with B1map cases */
t1w_post.into{t1w_post_a;t1w_post_b}

/* Split mtsat_from_preproc into two to deal with B1map cases */
mtsat_from_preproc.into{mfp_a;mfp_b}

/*Merge tw1_post with mtsat_from_preproc and b1plus.*/
t1w_post_a
    .join(mfp_a)
    .join(b1plus)
    .set{mtsat_for_fitting_with_b1}

/*Merge tw1_post with mtsat_from_preproc only.*/
t1w_post_b
    .join(mfp_b)
    .set{mtsat_for_fitting_without_b1}

/* When use_b1map set to true, process won't take those w/o a matching B1map*/
process Fit_MTsat_With_B1map{
    cpus 2

    publishDir = "$root/derivatives/qMRLab"

    when:
        params.use_b1map == true

    input:
        set sid, file(t1w), file(mtw_reg), file(pdw_reg),\
        file(mask), file(b1map) from mtsat_for_fitting_with_b1


    output:
        file "${sid}_dene.txt"

    script: 
        """
        echo "$mtw_reg\n" >> ${sid}_dene.txt 
        echo "$pdw_reg\n" >> ${sid}_dene.txt
        echo "$t1w\n" >> ${sid}_dene.txt 
        echo "$mask\n" >> ${sid}_dene.txt 
        echo "$b1map\n" >> ${sid}_dene.txt 
        echo "With B1map\n" >> ${sid}_dene.txt  
        """
}

/*Agnostic to any B1maps that may exist if params.use_b1map set to false..*/
process Fit_MTsat_Without_B1map{
    cpus 2

    publishDir = "$root/derivatives/qMRLab"
    
    when:
        params.use_b1map == false

    input:
        set sid, file(t1w), file(mtw_reg), file(pdw_reg),\
        file(mask) from mtsat_for_fitting_without_b1

    output:
        file "${sid}_dene.txt"

    script: 
        """
        echo "$mtw_reg\n" >> ${sid}_dene.txt 
        echo "$pdw_reg\n" >> ${sid}_dene.txt
        echo "$t1w\n" >> ${sid}_dene.txt 
        echo "$mask\n" >> ${sid}_dene.txt 
        echo "Without B1map\n" >> ${sid}_dene.txt  
        """
}

