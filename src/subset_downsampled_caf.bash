#!/bin/bash
gzcat data/caf/full_lof_metrics_by_transcript_an_adj.downsamplings.txt.bgz | head -1 > data/caf/downsampled.subset.txt
gzcat data/caf/full_lof_metrics_by_transcript_an_adj.downsamplings.txt.bgz | grep global | grep 10824 >> data/caf/downsampled.subset.txt
