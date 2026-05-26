include { SEQKIT_SORT } from '../../../modules/nf-core/seqkit/sort/main'

workflow SEQKIT_FASTQ_SORT {

    take:
    ch_fastq

    main:

    // Sort fastq reads to make bwa deterministic
    // See https://www.biostars.org/p/238628/#238817
    // and https://github.com/malariagen/pipelines/issues/41
    // For bwa mem, there are two sources of randomness:
    // 1. insert size estimation, which depends on chunk size, which depends on number of threads (chunk size = n_threads * 10M)
    //  Using `-K 100000000` allows for deterministic results regardless of the number of threads.
    //  Supposedly also sensitive to order of reads?
    //  https://github.com/kaist-ina/BWA-MEME/issues/27
    // 2. for alignments with multiple alternate positions with the same score
    //  (i.e. those with XA tags), the primary alignment is assigned randomly
    //  and this depends on the order of the fastq reads.
    //  https://github.com/lh3/bwa/issues/192

    // Split paired-end reads into separate channel elements with their own metamap
    ch_fastq_unpaired = ch_fastq.flatMap { meta, reads ->

        // Single-end files: pass through untouched
        if ( meta.single_end ) {
            return [ tuple(meta, reads) ]
        }

        // Paired-end files: emit one tuple for each read pair [meta, fq_R#]
        reads.withIndex(1).collect { fq, read_num ->

                // Assign meta.id that includes read num to make seqkit_sort produce distinguishable file names
                tuple(
                    meta + [
                        pair_id: meta.id,
                        id: "${meta.id}_R${read_num}",
                        read_num: read_num
                    ],
                    fq)
            }
    }                                                                   // [[id:s_R1, pair_id:s-LANE_L00#-FC_####, single_end:false, sample:s, ..., read_num:1], s_R1.fq.gz]

    // Sort each fq file independently
    SEQKIT_SORT(ch_fastq_unpaired)

    // Restore original read-pair level id and reconstruct PE structure
    ch_fastq_sorted = SEQKIT_SORT.out.fastx                             // [[id:s_R1, pair_id:s-LANE_L00#-FC_####, single_end:false, sample:s, ..., read_num:1], s_R1.fastq.gz]
        .map { meta, fq ->

            def restored_meta = meta.pair_id ?
                meta + [ id: meta.pair_id ] :
                meta

            tuple(restored_meta, fq)
        }                                                               // [[id:s-LANE_L00#-FC_####, ...], s_R1.fastq.gz]
        .branch { meta, _fq ->
            pe: !meta.single_end
            se:  meta.single_end
        }

    ch_fastq_pe = ch_fastq_sorted.pe
        .map { meta, fq ->
            tuple(meta.id, [meta.read_num, fq], meta)
        }                                                               // [s-LANE_L00#-FC_####, [1, s_R1.fastq.gz], [meta]]
        .groupTuple()                                                   // [s-LANE_L00#-FC_####, [[1, s_R1.fastq.gz], [2, s_R2.fastq.gz]], [[meta_R1], [meta_R2]]
        .map { _key, read_pairs, metas ->

            def ordered_fqs = read_pairs
                .sort { it[0] }
                .collect { it[1] }

            def meta = metas[0].subMap(
                metas[0].keySet() - ['read_num', 'pair_id']
            )

            tuple(meta, ordered_fqs)
        }                                                               // [[id, single_end, num_entries, multiple_lanes, multiple_libraries, multiple_flowcells, read_group], [s_R1.fastq.gz, s_R2.fastq.gz]]

    // Final combined channel
    ch_fastq = ch_fastq_pe.mix(ch_fastq_sorted.se)

    emit:
    fastq_sorted = ch_fastq
}
