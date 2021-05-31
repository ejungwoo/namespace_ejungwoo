PWDESC=$(echo $PWD | sed 's_/_\\/_g')
sed "s/PWD/$PWDESC/" input/rootlogon.C > rootlogon.C 

echo ================================
echo - rootlogon configuration
echo   Copy :
echo '    ' gROOT -\> LoadMacro\(\"$PWD/rootlogon.C\"\)\;
echo   Paste to :
echo '    ' $KEBIPATH/macros/rootlogon.C
echo ================================
#echo - List of conf names :
#ls input/*.conf | xargs -n 1 basename | sed "s/.conf//"
#echo ================================
