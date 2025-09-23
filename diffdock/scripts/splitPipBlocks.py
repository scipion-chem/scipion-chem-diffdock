import os, sys

def parsePipDependencies(ymlFile):
  headLines, pipDic = '', {}
  inPipDeps, pipDepIdx = False, 0

  with open(ymlFile) as f:
    for line in f:
      if line.strip() == '- pip:':
        inPipDeps = True
        pipDepIdx += 1

      if inPipDeps:
        if not pipDepIdx in pipDic:
          pipDic[pipDepIdx] = ''
        pipDic[pipDepIdx] += line
      else:
        headLines += line
  return headLines, pipDic


def splitPipDependencies(ymlFile):
  '''Method to break down the different pip installation blocks DiffDock (wrongly) offers in their YML
  Returns the list of yml files that must be installed sequentially (first using "conda env create -f file1.yml "
  and the rest "conda env update -f filex.yml"
  '''
  headLines, pipDic = parsePipDependencies(ymlFile)

  envFiles = []
  envFileBase = os.path.join(os.path.dirname(ymlFile), 'env_{}.yml')
  for idx, pipDeps in pipDic.items():
    envFiles.append(envFileBase.format(idx))
    with open(envFiles[-1], 'w') as fo:
      fo.write(headLines)
      fo.write(pipDeps)
  return envFiles

if __name__ == "__main__":
  ymlFile = sys.argv[1]
  retFirst = True if (len(sys.argv) > 2 and sys.argv[2].lower()) == 'true' else False
  envFiles = splitPipDependencies(ymlFile)
  if retFirst:
    print(envFiles[0])
  else:
    for f in envFiles[1:]:
      print(f)