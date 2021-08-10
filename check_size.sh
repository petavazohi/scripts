#/bin/bash
du -s * | sort -rn | cut -f2- | xargs -d "\n" du -sh
