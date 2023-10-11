# Improving-EHBSA-with-Relational-Decomposition-and-Entropy-aware-Sampling

thesis: https://drive.google.com/file/d/1Duz9avtO2I9WyoUE-s92sSn6wuRzB0qH/view?usp=sharing

oral presentation: https://drive.google.com/file/d/1YZcnGlPHuk0G3v1-t_MEoJY9v65x5nKC/view?usp=sharing

### usage
```
$ make
$ ./RankingEDAsCEC -i [problem_instance_file_name]
                   -o [result_file_name]
                   -s [random_seed]
                   -t [problem_type]
                   -m [model_name]
                   -d C -v 0
                   -p [population_size]
                   -e [max_NFE]
                   -b [bias_ratio]
                   -c [number_of_cut_point]
```

example:
```
$ make
$ ./RankingEDAsCEC -i TSP_instance/gr24.tsp
                   -o result.txt
                   -s 0
                   -t TSP
                   -m NRBOP
                   -d C -v 0
                   -p 200
                   -e 57600
                   -b 0.0002
                   -c 3
```

installed problem: TSP, QAP, PFSP, LOP, SOP
