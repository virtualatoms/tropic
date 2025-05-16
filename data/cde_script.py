from os import listdir
from chemdataextractor import Document
from chemdataextractor.model.model import Compound


if __name__ == "__main__":
    papers = listdir("./literature")
    for paper in papers:
        with open(f"./literature/{paper}", "rb") as fstream:

            # use HTML or XML instead?
            doc = Document.from_file(fstream)
            doc.models = [Compound]
            for record in doc.records:
                print(record.serialize())
