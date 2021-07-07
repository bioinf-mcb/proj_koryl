import os

def folder(parent_dir, folder_name):
    '''
    This function create new folder in directory
    '''
    path = os.path.join(parent_dir, folder_name)
    os.mkdir(path)
    print("Directory '% s' created" % folder_name)
    
