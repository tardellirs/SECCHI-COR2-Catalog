pro clickjmap_cor2_test,jmap,spfile=spfile

  ;+
  ;
  ;   Name: clickjmap_cor2_test
  ;
  ;   Purpose: Interface for select points in the HTmap (Clicks)
  ;
  ;-

  x=fltarr(2)
  y=fltarr(2)
  xp=x & yp=y
  wid=!D.Window                                       ; id = "current window"
  WSet, wid
  xsize = float(!D.X_VSize)                           ; X and Y size of the current visible area
  ysize = float(!D.Y_VSize)                           ;
  Window, /Pixmap, /Free, XSize=xsize, YSize=ysize    ; Set a window for displaying.
  pixID = !D.Window
  Device, Copy=[0, 0, xsize, ysize, 0, 0, wid]        ; Copy window contents
  WSet, wid                                           ; Work in the new window, Here we can draw
  wshow,wid                                           ; on the figure without changing the data.
  
  i=0
  a=1
  while a eq 1 do begin
    cursor,xvar,yvar,/change                        ; Get any cursor movement
    !mouse.button=0                                   ; Clear mouse button flag
    cursor,xvar,yvar,/change                      ; Check any cursor movement
    if !mouse.button EQ 1 THEN BEGIN                  ; If the left mouse button was pressed,
      i++
      xjmap = xvar                                 ; Concatenate X00 array
      yjmap = yvar                                ; Concatenate Y00 array
      acc_xp = exist(acc_xp) ? [acc_xp,xjmap] : xjmap
      acc_yp = exist(acc_yp) ? [acc_yp,yjmap] : yjmap
      xp[i-1]=xjmap
      yp[i-1]=yjmap
      x[i-1]=xjmap
      y[i-1]=yjmap
      wait,0.1
    endif
    if exist(acc_xp) && !mouse.button EQ 0 then begin
      acc_xp_tmp = [acc_xp,xvar]
      acc_yp_tmp = [acc_yp,yvar]
    endif else begin
      acc_xp_tmp = xvar
      acc_yp_tmp = yvar
    endelse
    Device, Copy=[ 0, 0, xsize, ysize, 0, 0, pixID]   ; Clear the line drawn.
    plots,acc_xp_tmp,acc_yp_tmp,psym=6,SYMSIZE=.5,color=255,thick=1
    if i eq 2 then begin
      i=0
      for j=0,n_elements(acc_xp)-2,2 do plots,acc_xp[j:j+1],acc_yp[j:j+1],color=255 ;; aa don't exist.
      if not exist(speedfile) then begin
        speedfile=float([x,y])       ; Format: X0  X1  Y0  Y1
      endif else begin
        speedfile=[[speedfile],[x,y]]
      endelse
      ind_xjmap = exist(ind_xjmap)?[ind_xjmap,xjmap]:xjmap
      ind_yjmap = exist(ind_yjmap)?[ind_yjmap,yjmap]:yjmap
      ind_numb = exist(ind_numb)?[ind_numb,n_elements(acc_xp)/2]:(n_elements(acc_xp)/2)
    endif else begin
      for j=0,n_elements(acc_xp_tmp)-2,2 do plots,acc_xp_tmp[j:j+1],acc_yp_tmp[j:j+1],color=255
    endelse
    if !mouse.button eq 4 then a=2
  endwhile
  
  if exist(acc_xp) then begin
    if i ne 0 then begin
      acc_xp=[acc_xp,xvar] & acc_yp=[acc_yp,yvar]
      xj=fltarr(2) & yj=fltarr(2) & x=fltarr(2) & y=fltarr(2)
      xj[1]=acc_xp[n_elements(acc_xp)-1] & xj[0]=acc_xp[n_elements(acc_xp)-2]
      yj[1]=acc_yp[n_elements(acc_yp)-1] & yj[0]=acc_yp[n_elements(acc_yp)-2]
      x[1]=round(xj[1])
      y[1]=round(yj[1])
      x[0]=round(xj[0])
      y[0]=round(yj[0])
      speedfile = exist(speedfile) ? [[speedfile],[x,y]] : float([x,y])
    endif
    Device, Copy=[ 0, 0, xsize, ysize, 0, 0, pixID]
    plots,acc_xp,acc_yp,psym=6,SYMSIZE=.5,color=255,thick=1
    for j=0,n_elements(acc_xp)-2,2 do plots,acc_xp[j:j+1],acc_yp[j:j+1],color=255
  endif
  spfile=speedfile
end