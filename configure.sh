PWDESC=$(echo $PWD | sed 's_/_\\/_g')

if [[ -v KEBIPATH ]]
then
  sed -e "s/PWD/$PWDESC/" -e "s/KEBIPATHISSET/0/" input/rootlogon.C > rootlogon.C 
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
else
  sed -e "s/PWD/$PWDESC/" -e "s/KEBIPATHISSET/1/" input/rootlogon.C > rootlogon.C 
  echo ================================
  echo - rootlogon configuration
  echo '  echo' Rint.Logon: $PWD\/rootlogon.C \>\> ~\/.rootrc
  echo ================================
fi
