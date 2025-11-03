.PHONY: all analysis plots report quick-demo montecarlo latex clean

all: analysis plots report

analysis:
	Rscript scripts/comprehensive_study.R

plots:
	Rscript scripts/create_visualizations.R

report:
	Rscript scripts/generate_report.R

quick-demo:
	Rscript scripts/arima_oil_quick_demo.R

montecarlo:
	Rscript scripts/run_monte_carlo.R

latex:
	latexmk -pdf -cd latex/PMM2_ARIMA.tex

clean:
	find results -mindepth 1 ! -name '.gitignore' -delete
