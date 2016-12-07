echo Building Analyzer
if [ ! -d "asm" ]; then
  echo Configuring Analyzer
  mkdir asm
  cd asm && emconfigure ../configure --disable-multithread --disable-runtime-cpu-detect --target=generic-gnu --enable-accounting
fi

cd asm
emmake make
cp examples/analyzer_decoder examples/analyzer_decoder.bc
emcc -O3 examples/analyzer_decoder.bc -o examples/decoder.js -s TOTAL_MEMORY=134217728 -s MODULARIZE=1 -s EXPORT_NAME="'DecoderModule'" --post-js "../ins/post.js" --memory-init-file 0
cd ..
mkdir -p ins/bin
cp asm/examples/decoder.js ins/bin/decoder.js
tsc
echo Analyzer is ready, serve it from the ins directory using your favorite web server. E.g. 'cd ins && python -m SimpleHTTPServer'
