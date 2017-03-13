
# An example how to use snakemake's wildcards


#just create example files for input
open("B.se.fq", "w").close()
open("A.pe1.fq", "w").close()	
open("A.pe2.fq", "w").close()	
#------


configfile: "config.json"


rule all:
    input:
        expand("{dataset}.sorted.txt",dataset=config["samples"])


rule sort:
    input:
        lambda wildcards: config["units"][wildcards.harry]
    output:
        "{harry}.sorted.txt"
    log:
        "logs/sort/{harry}.log"
    shell:
        "sort {input} > {output};"
        "cat {input};"
        "echo {wildcards.harry};"
        "echo '--------';"
        "cat {output};"
        "2> {log}"
