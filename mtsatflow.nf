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

Users: Please see USAGE for usage details.
 */

/*Set defaults for parameters determining logic flow to false*/
params.root = false 
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

if(params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    
    /* ==== Custom: MTSat inputs ==== */  
    /* Here, alphabetical indexes matter. Therefore, MTw -> PDw -> T1w */
    in_data = Channel.fromFilePairs("$root/**/*{MTw.nii.gz,PDw.nii.gz,T1w.nii.gz,mt_sat_prot.json}", maxDepth: 1, size: 4, flat: true){it.parent.name}

    (mtw, pdw, t1w, json) = in_data
        .map{sid, MTw, PDw, T1w, mt_sat_prot  -> [tuple(sid,MTw),      
                                         tuple(sid, PDw),
                                         tuple(sid, T1w),
                                         tuple(sid, mt_sat_prot)]}                                   
        .separate(4)
    
    /* ==== Custom: B1 map ==== */  
    /* Look for B1map in subject folders */
    Channel 
    .fromPath("$root/**/*B1plusmap.nii.gz",
                    maxDepth:1)                              
    .map{[it.parent.name, it]}
    .set{b1plus}

}
else{
    error "ERROR: Argument (--root) must be passed. See USAGE."
}

/*Each data type is defined as a channel. To pass all the channels 
  to the same process accurately, these channels must be joined. 
*/ 

/*Split T1w into three channels
    t1w_pre_ch1 --> mtsat_for_alignment
    t1w_pre_ch2 --> mtsat_for_bet
    t1w_pre_ch3 --> t1w_post
*/
t1w.into{t1w_pre_ch1; mtsat_for_bet; t1w_post}

/*After alignment, we still need simpleName of the parent files 
 access .json content. 
    - Given the B1map and Mask options, there are 4 possible processes.
      Original mtw, pdw and t1w channels should be accessible by them.    
*/


/* Merge PDw, MTw and T1w for alignment*/
pdw 
    .join(mtw)
    .join(t1w_pre_ch1)
    .set{mtsat_for_alignment}

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
log.info "Dimensionality: $params.ants_dim"
log.info "Metric: $params.ants_metric"
log.info "Weight: $params.ants_metric_weight"
log.info "Number of bins: $params.ants_metric_bins"
log.info "Sampling type: $params.ants_metric_sampling"
log.info "Sampling percentage: $params.ants_metric_samplingprct"
log.info "Transform: $params.ants_transform"
log.info "Convergence: $params.ants_convergence"
log.info "Shrink factors: $params.ants_shrink"
log.info "Smoothing sigmas: $params.ants_smoothing"
log.info ""
log.info "[qMRLab mt_sat]"
log.info "---------------"
log.info ""

if (params.USE_B1){
log.info "B1+ correction has been ENABLED."  
log.warn "Process will be skipped for participants missing a B1map file."   
log.info "B1 correction factor: $params.COR_B1"}
if (!params.USE_B1){
log.info "B1+ correction has been DISABLED."
log.warn "Process will NOT take any (possibly) existing B1maps into account."
}

/*Perform rigid registration to correct for head movement across scans:
    - MTw (moving) --> T1w (fixed)
    - PDw (moving) --> T1w (fixed)
*/     

process Align_Input_Volumes {

    publishDir = "$root/derivatives/qMRLab/" 

    input:
    set sid, file(pdw), file(mtw), file(t1w) from mtsat_for_alignment

    output:
    set sid, "${sid}_acq-MTon_MTS_aligned.nii.gz", "${sid}_acq-MToff_MTS_aligned.nii.gz"\
    into mtsat_from_alignment
   
    script:
    """
    touch ${sid}_acq-MTon_MTS_aligned.nii.gz
    touch ${sid}_acq-MToff_MTS_aligned.nii.gz
    """
}

process Extract_Brain{

    publishDir = "$root/derivatives/qMRLab/"

    when:
        params.USE_BET == true

    input:
    set sid, file(t1w) from mtsat_for_bet

    output:
    set sid, "${sid}_acq-T1w_mask.nii.gz" optional true into mtsat_from_bet
    
    script:
         if (params.bet_recursive){
            """    
            touch ${sid}_acq-T1w_mask.nii.gz
            """}
        else{
            """    
            touch ${sid}_acq-T1w_mask.nii.gz
            """
        }

}

/* Split t1w_post into two to deal with B1map cases */
t1w_post.into{t1w_post_ch1;t1w_post_ch2}

/* Split mtsat_from_alignment into two to deal with B1map cases */
mtsat_from_alignment
    .join(json)
    .into{mfa_ch1;mfa_ch2}

/* There is no input optional true concept in nextflow
The process consuming the individual input channels will 
only execute if the channel is populated.
*/

/* We need to conditionally create channels with
input data or as empty channels.
*/
if (!params.USE_BET){
    Channel
        .empty()
        .set{mtsat_from_bet}
}

/* Split mtsat_from_bet into two to deal with B1map cases later. */
mtsat_from_bet.into{mtsat_from_bet_ch1;mtsat_from_bet_ch2}

/*Merge tw1_post with mtsat_from_alignment and b1plus.*/
t1w_post_ch1
    .join(mfa_ch1)
    .join(b1plus)
    .set{mtsat_for_fitting_with_b1}

mtsat_for_fitting_with_b1.into{mtsat_with_b1_bet;mtsat_with_b1}

/*Merge tw1_post with mtsat_from_alignment only.
WITHOUT B1 MAP
*/
t1w_post_ch2
    .join(mfa_ch2)
    .set{mtsat_for_fitting_without_b1}

mtsat_for_fitting_without_b1.into{mtsat_without_b1_bet;mtsat_without_b1}

/* We need to join these channels to avoid problems.
WITH B1 MAP
*/
mtsat_with_b1_bet
    .join(mtsat_from_bet_ch1 )
    .set{mtsat_with_b1_bet_merged}

/* Depeding on the nextflow.config 
settings for USE_B1 and USE_BET, one of th
following 4 processes will be executed. 
*/

process Fit_MTsat_With_B1map_With_Bet{

    publishDir = "$root/derivatives/qMRLab/"

    when:
        params.USE_B1 == true && params.USE_BET == true

    input:
        set sid, file(t1w), file(mtw_reg), file(pdw_reg), file(json),\
        file(b1map), file(mask) from mtsat_with_b1_bet_merged
        

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_MTsat.nii.gz" 
        file "${sid}_mt_sat.qmrlab.mat"

    script: 
        if (params.PLATFORM == 'octave'){
                log.info "qMRLab::mt_sat | Octave"
                """
                    octave --no-gui --eval "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'mask','$mask','b1map','$b1map','b1factor',$params.COR_B1,'custom_json','$json')"
                """
                } else{
                log.info "qMRLab::mt_sat | MATLAB"    
                """
                  matlab -nodisplay -nosplash -nodesktop -r "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'mask','$mask','b1map','$b1map','b1factor',$params.COR_B1,'custom_json','$json')"
                """
                }
            
}

process Fit_MTsat_With_B1map_Without_Bet{
    cpus 1

    publishDir = "$root/derivatives/qMRLab/"

    when:
        params.USE_B1 == true && params.USE_BET == false

    input:
        set sid, file(t1w), file(mtw_reg), file(pdw_reg), file(json),\
        file(b1map) from mtsat_with_b1

    output:
       file "${sid}_T1map.nii.gz"
       file  "${sid}_MTsat.nii.gz"
       file "${sid}_mt_sat.qmrlab.mat"

    script: 
             if (params.PLATFORM == 'octave'){
                log.info "qMRLab::mt_sat | Octave"
                """
                    octave --no-gui --eval "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'b1map','$b1map','b1factor',$params.COR_B1,'custom_json','$json')"
                """
                } else{
                log.info "qMRLab::mt_sat | MATLAB"    
                """
                  matlab -nodisplay -nosplash -nodesktop -r "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'b1map','$b1map','b1factor',$params.COR_B1,'custom_json','$json')"
                """
                }
            
}


/* We need to join these channels to avoid problems.*/
mtsat_without_b1_bet
    .join(mtsat_from_bet_ch2)
    .set{mtsat_without_b1_bet_merged}

process Fit_MTsat_Without_B1map_With_Bet{
    cpus 1

    publishDir = "$root/derivatives/qMRLab/"
    
    when:
        params.USE_B1 == false && params.USE_BET==true

    input:
        set sid, file(t1w), file(mtw_reg), file(pdw_reg), file(json),\
        file(mask) from mtsat_without_b1_bet_merged

    output:
        file "${sid}_T1map.nii.gz"
        file "${sid}_MTsat.nii.gz"
        file  "${sid}_mt_sat.qmrlab.mat"

    script: 
             if (params.PLATFORM == 'octave'){
                log.info "qMRLab::mt_sat | Octave"
                """
                    octave --no-gui --eval "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'mask','$mask','custom_json','$json')"
                """
                } else{
                log.info "qMRLab::mt_sat | MATLAB"
                """
                  matlab -nodisplay -nosplash -nodesktop -r "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'mask','$mask','custom_json','$json')"
                """
                }
            
}

process Fit_MTsat_Without_B1map_Without_Bet{
    cpus 1

    publishDir = "$root/derivatives/qMRLab/"
    
    when:
        params.USE_B1 == false && params.USE_BET==false

    input:
        set sid, file(t1w), file(mtw_reg), file(pdw_reg),\
        file(json) from mtsat_without_b1

    output:
        file "${sid}_T1map.nii.gz"
        file  "${sid}_MTsat.nii.gz"
        file  "${sid}_mt_sat.qmrlab.mat"

    script: 
            if (params.PLATFORM == 'octave'){
                log.info "qMRLab::mt_sat | Octave"
                """
                    octave --no-gui --eval "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'custom_json','$json')"
                """
                } else{
                log.info "qMRLab::mt_sat | MATLAB"    
                """
                  matlab -nodisplay -nosplash -nodesktop -r "mt_sat_wrapper('${mtw_reg.simpleName}.nii.gz','${pdw_reg.simpleName}.nii.gz','${t1w.simpleName}.nii.gz',[],[],[],'custom_json','$json')"
                """
                }
            
}

