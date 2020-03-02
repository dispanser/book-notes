SUB=$1
# ls -l $SUB/[09][09]*.md
pandoc $SUB/[09][09]*.md --highlight-style=tango -s -o "$SUB.pdf"
