def get_sub_list(DATA_DIR):
    # Get the list of files in the current directory
    files = os.listdir(DATA_DIR)
    
    # Sort the files by modification time
    sorted_files = sorted(files)
    
    sub_list = []

    for i, filename in enumerate(sorted_files):
        if filename[0:3] == 'sub':
            sub = filename
            sub_folder = osp.join(DATA_DIR, filename)
            ses_count = 0
            ses_nums = []
            
            for sub_filename in os.listdir(sub_folder):
                if sub_filename[0:3] == 'ses':
                    ses_count += 1
                    ses_nums.append(sub_filename)
                    
            for s in range(ses_count):
                sub_list.append((sub,ses_nums[s]))