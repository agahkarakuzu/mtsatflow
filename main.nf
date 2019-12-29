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
    .fromPath("$root/**/fmap/sub-*_B1plusmap.nii.gz",
                    maxDepth:2)                
    .map{[it.parent.parent.name, it]}
    .into{b1plus}   
}   
/*If data is not BIDS-compatible, define a custom directory structure
    Please see USAGE for further details. 
*/ 
else if(params.root && !params.bids){

}
else{
    error "ERROR: Arguments (--root) and (--bids) must be passed. See USAGE."
}

/*Each data type is defined as a channel. To pass all the channels 
  to the same process accurately, these channels must be joined. 
*/ 
pdw.set{pdw_ch}

pdw_ch
    .join(mtw)
    .join(t1w)
    .set{mtsat_ch}

log.info "qMRflow: MTsat pipeline"
log.info "======================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.info ""
log.info "DATA"
log.info ""
if (params.bids){
log.info "== BIDS option has been enabled."
log.info "qMRI protocols will be read from sidecar .json files."
log.info "If a B1plusmap is provided in the fmap folder, correction factor defined in the nexflow.config will be taken into account."
}
else{
log.info "== Custom file/folder organization described in USAGE will be assumed."
log.info "If an mtsat_protocol.json file is provided for a subject, acquisition metadata will be overridden. Otherwise, those defined in nextflow.config will be used."
}

log.info ""
log.info "OPTIONS"
log.info "======="
log.info ""
log.info "[ANTs Registration]"
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
log.info ""
log.info "NOTE: mt_sat protocol values are subjected to change based on the (--bids) option. See more at USAGE."
log.info ""
log.info "Default values provided in the nextflow.config: "
log.info ""
log.info "Flip angles:\n\t MTon: $params.mtw_fa\n\t MToff: $params.pdw_fa\n\t T1w: $params.t1w_fa"
log.info ""
log.info "Repetition times:\n\t MTon: $params.mtw_tr\n\t MToff: $params.pdw_tr\n\t T1w: $params.t1w_tr"
log.info ""
log.info "B1 correction factor: $params.b1_cor_factor"


/*Perform rigid registration to correct for head movement across scans:
    - MTw (moving) --> T1w (fixed)
    - PDw (moving) --> T1w (fixed)
*/     
/*   
process Align_Input_Volumes {
    cpus 2
    container 'qmrlab/ants'

    if(params.root && params.bids){
        publishDir = "$root/derivatives/qMRLab"
    }

    input:
    set sid, file(pdw), file(mtw), file(t1w) from mtsat_ch

    output:
    set sid, "${sid}_acq-MTon-to-T1w_MTS.nii.gz", "${sid}_acq-MToff-to-T1w_MTS.nii.gz"\
        into mtsat_ch_reg

    script:
    """
    antsRegistration -d $params.dim\\ 
        --float 0\\ 
        -o [${sid}_mtw2t1w,${sid}_acq-MTon-to-T1w_MTS.nii.gz]\\ 
        --transform $params.transform\\ 
        --metric $params.metric[$t1w,$mtw,$params.metric_weight,\\
        $params.metric_bins,$params.metric_sampling,$params.metric_samplingprct]\\ 
        --convergence $params.convergence\\ 
        --shrink-factors $params.shrink\\ 
        --smoothing-sigmas $params.smoothing

    antsRegistration -d $params.dim\\ 
        --float 0\\ 
        -o [${sid}_pdw2t1w,${sid}_acq-MToff-to-T1w_MTS.nii.gz]\\ 
        --transform $params.transform\\ 
        --metric $params.metric[$t1w,$pdw,$params.metric_weight,\\
        $params.metric_bins,$params.metric_sampling,$params.metric_samplingprct]\\ 
        --convergence $params.convergence\\ 
        --shrink-factors $params.shrink\\ 
        --smoothing-sigmas $params.smoothing

    """
}
*/    