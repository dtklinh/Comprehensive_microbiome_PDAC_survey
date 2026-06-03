rule: Chap1_Alpha
    input: "data/pseq.rds"
    output: 
        "results/Chap1/.png"
    conda: "envs/Chap1.yaml"
    shell: """
        Rscript scripts/Chap1.B.py {input}
    """