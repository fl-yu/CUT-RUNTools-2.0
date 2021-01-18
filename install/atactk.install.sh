#the current directory should contain make_cut_matrix.patch, and metrics.py.patch
cur=`pwd`

mkdir git
cd git
git clone https://github.com/ParkerLab/atactk
cd atactk
git checkout 6cd7de0
cd ..
pip install --user ./atactk
cd ..

cd ~/.local/bin
cp $cur/make_cut_matrix.patch .
patch -p0 -N --dry-run --silent make_cut_matrix < make_cut_matrix.patch 2> /dev/null
if [ $? -eq 0 ];
then
patch -p0 -N make_cut_matrix < make_cut_matrix.patch
fi

cd ~/.local/lib/python2.7/site-packages/atactk
cp $cur/metrics.py.patch .
patch -p0 -N --dry-run --silent metrics.py < metrics.py.patch 2> /dev/null
if [ $? -eq 0 ];
then
patch -p0 -N metrics.py < metrics.py.patch
fi

cd $cur
