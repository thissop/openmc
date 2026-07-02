import os

MB_THRESHOLD = 75

# Size threshold in bytes (MB_THRESHOLD MB)
SIZE_THRESHOLD = MB_THRESHOLD * 1024 * 1024
gitignore_path = '.gitignore'

large_files = []

# Walk through the repository directory
for root, dirs, files in os.walk('.'):
    # Skip the .git directory to avoid modifying Git's internal state
    if '.git' in dirs:
        dirs.remove('.git')
        
    for file in files:
        file_path = os.path.join(root, file)
        
        # Get file size
        try:
            file_size = os.path.getsize(file_path)
        except OSError:
            continue  # Skip broken symlinks or files with permission issues
            
        if file_size > SIZE_THRESHOLD:
            # Get relative path and standardize separators for .gitignore (forward slashes)
            rel_path = os.path.relpath(file_path).replace(os.path.sep, '/')
            large_files.append(rel_path)

# Write to .gitignore
if large_files:
    # Read existing entries to prevent duplicates
    existing_entries = set()
    if os.path.exists(gitignore_path):
        with open(gitignore_path, 'r') as f:
            existing_entries = {line.strip() for line in f if line.strip()}
            
    with open(gitignore_path, 'a') as f:
        # Add a header if the file is new
        if not existing_entries:
            f.write(f"# Automatically ignored large files (> {MB_THRESHOLD} MB)\n")
            
        count = 0
        for path in large_files:
            if path not in existing_entries:
                f.write(f"/{path}\n")
                count += 1
                
        print(f"Successfully added {count} files larger than {MB_THRESHOLD} MB to .gitignore.")
else:
    print(f"No files larger than {MB_THRESHOLD} MB found.")