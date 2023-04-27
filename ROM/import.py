from foamfile import FoamFile

with FoamFile("") as f:
    foam_content = f.read()
    print(f.header)
    print(foam_content)

with FoamFile("path/to/file", "w", foam_class="dictionary") as f:
    f.write(foam_content)