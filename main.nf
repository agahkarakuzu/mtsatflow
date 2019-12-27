#!/usr/bin/env nextflow

params.root = false 
params.bids = false 
params.help = false

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


if(params.root && params.bids){
    log.info "Input: $params.root"
    root = file(params.root)
    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_acq-{MToff,MTon,T1w}_MTS.nii.gz", maxDepth: 2, size: 3, flat: true)

    (pdw, mtw, t1w) = in_data
        .map{sid, MToff, MTon, T1w  -> [    tuple(sid, MToff),
                                            tuple(sid, MTon),
                                            tuple(sid, T1w)]}                                   
        .separate(3)
       
}
else{
    error "pass --root"
}

pdw.set{pdw_ch}

pdw_ch
    .join(mtw)
    .join(t1w)
    .set{mtsat_ch}

log.info "MTsat pipeline"
log.info "==================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.info "Options"
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
        
process Align_Input_Volumes {
    cpus 2

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