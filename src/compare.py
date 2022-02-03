from filecmp import dircmp
def print_diff_files(dcmp):
  for name in dcmp.diff_files:
    print("diff_file %s found in %s and %s" % (name, dcmp.left,
              dcmp.right))
  for sub_dcmp in dcmp.subdirs.values():
    print_diff_files(sub_dcmp)
dcmp = dircmp('.', '/Users/heiland/git/pc4covid19-rheiland-master/src/custom_modules') 
print_diff_files(dcmp)
