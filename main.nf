params.root = false 
params.bids = false 

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


process Agah {
    cpus 1
    
     if(params.root && params.bids){
        publishDir = "$root/derivatives/qMRLab"
    }

    input:
    set sid, file(pdw), file(mtw), file(t1w) from mtsat_ch

    output:
    file "${sid}_readagah.txt"

    script:
    """
    echo "MToff --> ${sid} $pdw \n">> ${sid}_readagah.txt
    echo "MTon --> ${sid} $mtw \n">> ${sid}_readagah.txt
    echo "T1w --> ${sid} $t1w \n">> ${sid}_readagah.txt
    """
}