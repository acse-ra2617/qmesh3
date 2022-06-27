all: pdfs pngs

pdfs: architectureSchematic.pdf
pngs: architectureSchematic.png

%.pdf: %.svg
	inkscape $< -D -A $@ --export-text-to-path

%.png: %.svg
	inkscape $< -D --export-png=$@ --export-dpi=100

clean:
	@rm -f architectureSchematic.pdf architectureSchematic.png
