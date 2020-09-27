import os

os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("docs", exist_ok=True)

rule all:
	input: "docs/sim1.html"

rule data:
	input: "data/Living Circles Count Update.xlsx"
	output: "data/circles.rdata"
	shell:
		"cd scripts; Rscript data.r"

rule sim1:
	input: "data/circles.rdata"
	output: "results/sim1.rdata"
	shell:
		"cd scripts; Rscript sim1.r"

rule sim1_doc:
	input: "results/sim1.rdata", "docs/sim1.Rmd"
	output: "docs/sim1.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"sim1.Rmd\", output_format=\"all\")'"
