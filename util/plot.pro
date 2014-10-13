datadir='./BIN/'

for output=0,0 do begin
  ;fname=datadir+"Coldens."+string(output,format='(I4.4)')+".bin"
  fname="./BIN/test.bin"
  data=dblarr(256,256)
  openr,1,fname,/f77_unformatted
  readu,1,data
  close,1
  print,minmax(data)
  window,xsize=500,ysize=500
  loadct,33
  ;tvscale,alog10(data)
  tvscale,data
  ;set_plot,'ps'
  ;fname="Coldens."+string(output,format='(I4.4)')+".ps"
  ;device,file=fname,/color,bits=8,xsize=3,ysize=3,/encapsulated
  ;tvscale,alog10(data)
  ;device,/close
  ;set_plot,'x'
endfor


end
