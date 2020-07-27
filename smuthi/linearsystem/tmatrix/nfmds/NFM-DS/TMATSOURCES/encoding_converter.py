import io, glob

for src_path in glob.glob("*.f90"):
    dst_path='win/'+src_path
    print(dst_path)
    with io.open(src_path, mode="r", encoding="utf8") as fd:
        content = fd.read()
    with io.open(dst_path, mode="w", encoding="cp1252") as fd:
        fd.write(content)
