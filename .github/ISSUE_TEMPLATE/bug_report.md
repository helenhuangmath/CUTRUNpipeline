---
name: Bug report
about: A step failed or gave a wrong result
labels: bug
---

**What happened**
A short description of the problem.

**Command**
```
bash run_pipeline.sh config.sh ...
```

**Status and log**
Output of `bash check_status.sh config.sh`, and the last lines of the log with the error
(`Data/Log/*.err`):
```

```

**Settings**
- Pipeline version (`cat VERSION`):
- Cluster: LSF / SLURM / local
- Software: conda environments from `envs/` / modules
- Relevant lines of `config.sh` and the first lines of `sample_list.txt` (remove private paths):
```

```

**Expected result**
What you expected instead.
