# GitHub push notes

```bash
cd /Users/mark/github/STXit
unzip /path/to/STXit_final_merged_repo.zip
rsync -av STXit_final_merged_repo/ ./
rm -rf STXit_final_merged_repo
git add .
git commit -m "Replace scaffold with final merged STXit snapshot"
git push
```
