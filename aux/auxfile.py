def write4(path, arr1, arr2, arr3, arr4):
   
    target_file = open(path, 'w')

    for ind in range(len(arr1)):

        form = "{0:15.7e} {1:15.7e} {2:15.7e} {3:15.7e}\n".format(arr1[ind], \
                                                                  arr2[ind], \
                                                                  arr3[ind], \
                                                                  arr4[ind])

        target_file.write(form)

    target_file.close()

def num_lines(f):

    return sum(1 for line in open(f))
