### tc_qtclust*
script that runs QT clustering
requires Statistics::Basic module:
    http://search.cpan.org/~jettero/Statistics-Basic/

### ML_only.csv
expression values for Mnemiopsis leidyi genes during development

### colloblast_candidates.csv
list of colloblast candidate genes w corresponding ratio of expression before
and after tentacle time-window

### colloblast_candidates.clusters
clusters of colloblast candidates output of 
    `./tc_qtclust colloblast_candidates.csv > colloblast_candidates.clusters`

### tentacle_candidates.clusters
list of tentacle candidate genes w corresponding ratio of expression before
and after tentacle time-window

### tentacle_candidates.csv
clusters of tentacle candidates output of 
    `./tc_qtclust tentacle_candidates.csv > tentacle_candidates.clusters`


