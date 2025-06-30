import sys, os


def destrand( input_f, output_f):
   mw = open(output_f, 'w');
   
   with open( input_f, 'r') as mr:
       cache_line = [];
       while True:
           line = mr.readline();
           if not line: break;
           line = line.strip();
           l_lsp = line.split();
           
           if l_lsp[5]=='-':
              if len(cache_line)==0:
                  mw.write("{}\n".format("\t".join(l_lsp)) )
              else:
                  if l_lsp[0]==cache_line[0] and ( int(l_lsp[1])==int(cache_line[1])+1 ) and (cache_line[5]=='+' and l_lsp[5]=='-'):
                     cache_line[9] = str(int(l_lsp[9])+int(cache_line[9]) )
                     cache_line[10] = str((float(l_lsp[10])*float(l_lsp[9])  +float(cache_line[10])*float(cache_line[9]))/(int(l_lsp[9])+int(cache_line[9])) )
                     mw.write("{}\n".format("\t".join(cache_line)) )
                  else:
                     mw.write("{}\n".format("\t".join(cache_line)) )
                     mw.write("{}\n".format("\t".join(l_lsp)) )
                  cache_line = []
           elif  l_lsp[5]=='+':
              if not len(cache_line)==0:
                 mw.write("{}\n".format("\t".join(cache_line)) )
              cache_line = l_lsp   
           else:
              print("Warning!!! Unexpected: {}".format( l_lsp ))
   mw.close();



if __name__ == "__main__":
    input_file = sys.argv[1]
    ouput_file = os.path.dirname(input_file) + "/new_destranded_" + os.path.basename(input_file)
    if len(sys.argv) > 2: ouput_file = sys.argv[2]
    destrand( input_file, ouput_file)
    # ouput_file = sys.argv[2]
    # destrand( input_file, ouput_file)
    print(ouput_file, " done!")
