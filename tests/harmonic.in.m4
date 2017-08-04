harmonic {
    colvars        one
    centers        0.1
    forceConstant  0.001
ifdef(`centers_moving_full',include(centers-moving-full.in))`'ifdef(`centers_moving_half',include(centers-moving-half.in))`'ifdef(`centers_moving_stages',include(centers-moving-stages.in))}
