# Decoder error-path allocation hardening

`decode_frame_packed` previously retained an exception and traceback for every failed mask candidate on each corrupted strand. It now retains only the last failure, which preserves the raised `ValueError` and final exception cause while reducing allocations in noisy-mask fallback paths. This is a targeted code hygiene change, not a demonstrated reduction in 10 MiB indel timeouts. The existing test matrix and matched benchmark runs must validate behavioral equivalence before merging.
