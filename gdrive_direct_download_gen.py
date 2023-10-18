#!/usr/bin/env python3

import re
import argparse

def get_direct_download_link(link):
    pattern = r'id=([^&]+)'
    match = re.search(pattern, link)
    if match:
        file_id = match.group(1)
        # Construct the direct download link using the extracted file ID
        return f"https://drive.google.com/uc?export=download&id={file_id}"
    else:
        return None

if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Generate Google Drive direct download link')
    parser.add_argument('share_link', help='Google Drive share link')
    args = parser.parse_args()

    # Get the direct download link
    direct_download_link = get_direct_download_link(args.share_link)
    if direct_download_link:
        print(direct_download_link)
    else:
        print('Failed to generate direct download link')