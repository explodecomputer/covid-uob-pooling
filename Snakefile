import os

os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("docs", exist_ok=True)

rule all:
	input: "docs/sim.html"

rule data:
	input: "data/Living Circles Count Update.xlsx"
	output: "data/circles.rdata"
	shell:
		"cd scripts; Rscript data.r"

rule sim:
	input: "data/circles.rdata"
	output: "results/sim.rdata"
	shell:
		"cd scripts; Rscript sim.r"

rule sim_doc:
	input: "results/sim.rdata", "docs/sim.rmd"
	output: "docs/sim.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"sim.rmd\", output_format=\"all\")'"
