tail=_nonsyn
for filename in ./*variant_function; do
  id=`basename "${filename%.*}"`
  grep -E "nonsyn" $filename > tmp
  mv tmp $filename$tail 
done
