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

rule containment:
	input: "scripts/containment.r"
	output: "data/containment.rdata"
	shell:
		"cd scripts; Rscript containment.r"

rule ppv:
	input: "docs/ppv.rmd"
	output: "docs/ppv.html"

rule ct:
	input: "docs/ct.rmd"
	output: "data/efficiency_params.rdata", "docs/ct.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"ct.rmd\", output_format=\"all\")'"

rule lfd:
	input: "docs/lfd.Rmd"
	output: "data/lfd_fit.rdata", "docs/lfd.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"lfd.Rmd\", output_format=\"all\")'"

rule sim:
	input: "data/circles.rdata", "scripts/functions.r", "scripts/sim.r", "data/lfd_fit.rdata", "data/efficiency_params.rdata", "data/containment.rdata"
	output: "results/sim.rdata"
	shell:
		"cd scripts; Rscript sim.r"

rule sim_doc:
	input: "results/sim.rdata", "docs/sim.rmd"
	output: "docs/sim.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"sim.rmd\", output_format=\"all\")'"
