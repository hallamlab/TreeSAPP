#!/usr/bin/python

# MLTreeMap TSN v. 0.0
# FileHolder

class FileHolder:
    def __init__(self, _homeDirectory):
        self.listOfFiles = {}
        self.homeDirectory = str(_homeDirectory)

    def addFile(self, fileName, fileContents):
        self.listOfFiles{str(self.homeDirectory + fileName)} = str(fileContents)

    def readFile(self, fileName):
        if self.listOfFiles{str(self.homeDirectory + fileName)}:
            return fileContents
        else:
            raise error(str(self.homeDirectory + fileName) + ' not found!')

    def writeFile(self, fileName):
        try:
            out = open(str(self.homeDirectory + fileName), 'w')
            out.write(self.readFile(fileName))
            out.close()
        except:
            raise
