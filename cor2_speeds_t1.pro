pro COR2_speeds_T1, filename, inRadius, startTime=startTime, endTime=endTime, $
    xcen=xcen, ycen=ycen, angRange=angRange, sampleTime=sampleTime, runDiff=runDiff, $
    compare=compare
  ;+
  ;   Name: COR2_speeds_T1
  ;
  ;   Purpose: CMEs Speed determination for COR2 data
  ;
  ;   Input Parameters:
  ;      filename - MVI file name
  ;      inRadius - Radius of the circumference (Solar radius)
  ;
  ;   Keyword Parameters:
  ;      startTime - Start time  'yyyy-mm-dd hh:mm:ss'
  ;      endTime - End time  'yyyy-mm-dd hh:mm:ss'
  ;      xCen - Projection of the center on x-axis (pixel)
  ;      yCen - Projection of the center on y-axis (pixel)
  ;      angRange - Angle range [0.,360.]
  ;      runDiff - Running difference (number of images)
  ;      compare - Plot a comparison between Angle x Time and Distance x Time
  ;
  ;   Called routines:
  ;      ClickJmap_COR2_Test.pro
  ;      COR2Record_Test.pro
  ;      Refresh_Image.pro
  ;
  ;   Calling Sequence:
  ;      temporary9, filename, inRadius [,startTime=startTime, endTime=endTime, $
  ;                                    xCen=xCen, ycen=ycen, angRange=angRange, $
  ;                                    sampleTime=sampleTime, runDiff=runDiff ]
  ;
  ;   Calling Example:
  ;      COR2_speeds_T1, 'C:\sdodata\20090826\T_COR2_A_ratio.mvi',4
  ;      COR2_speeds_T1, 'C:\sdodata\20090826\T_COR2_A_ratio.mvi',4,sampleTime=300,/compare
  ;      
  ;      COR2_speeds_T1, 'C:\sdodata\20130101\test2.mvi',4
  ;
  ;      COR2_speeds_T1, 'C:\sdodata\20121001_05\test_cor2.mvi',4,startTime='2012-10-02 00:24:00',endTime='2012-10-04 00:24:00'
  ;      COR2_speeds_T1, '/Volumes/Disk1/WorkingOn/Teresa/MVIs/T_COR2_A_ratio.mvi', 4
  ;      COR2_speeds_T1, '/Volumes/Disk1/WorkingOn/Teresa/MVIs/T_COR2B_ratio.mvi', 4
  ;      COR2_speeds_T1, 'C:\sdodata\C2\T_C2_ratio.mvi',4
  ;      
  ;-
    
  ;----------------------------- Window size ----------------------------
    
  screen_resolution=get_screen_size()
  wsize_x = screen_resolution[0]*2./5.
  wsize_y = screen_resolution[1]*3./7.
  
  ;----------------- Read Header ---------------------------
  OPENR, lu, filename, /GET_LUN
  READ_MVI, lu, file_hdr, ihdrs, imgs, swapflag
  
  n_images = file_hdr.nf
  
  for i=0, n_images-1 do begin
    if i eq 0 then header = MVIHDR2STRUCT(ihdrs(i)) $
    else  header = [header, MVIHDR2STRUCT(ihdrs(i)) ]
  endfor
  
  detector = header[0].detector
  
  ;------------------ Get MVI filename and path -------------------------
  n_slash = n_elements(strsplit(filename,'\'))
  if n_slash eq 0 then begin
    n_slash = n_elements(strsplit(filename,'/'))
    mviFileName = strmid(filename,(strsplit(filename,'/'))[n_slash-1])
    mvipath = strmid(filename,0,(strsplit(filename,'/'))[n_slash-1])
  endif else begin
    mviFileName = strmid(filename,(strsplit(filename,'\'))[n_slash-1])
    mvipath = strmid(filename,0,(strsplit(filename,'\'))[n_slash-1])
  endelse
  
  ;------------------ Select the time range -----------------------------
  header_time = strarr(n_images)
  header_day = strarr(n_images)
  time_sec = fltarr(n_images)
  fileNameList=strarr(n_images)
  for i = 0,n_images-1 do begin
    header_time[i] = header[i].TIME_OBS
    header_day[i] = header[i].DATE_OBS
    fileNameList[i] = header[i].FILENAME
  endfor
  day = long(STRMID(header_day,8,2))
  hour = long(STRMID(header_time,0,2))
  min = long(STRMID(header_time,3,2))
  sec = long(round(float(STRMID(header_time,6,5))))
  time_sec = long(sec+min*60.+hour*60.*60.+day*60.*60.*24.)
  time_sec0 = time_sec
  if not exist(startTime) then begin
    startfile = 0 & endfile = n_images-1
    goto,noTimeInput
  endif
  startTimeSec = float(STRMID(startTime,8,2))*24.*60.*60.+float(STRMID(startTime,11,2))*60.*60.+ $
  float(STRMID(startTime,14,2))*60.+float(STRMID(startTime,17,2))
  endTimeSec = float(STRMID(endTime,8,2))*24.*60.*60.+float(STRMID(endTime,11,2))*60.*60.+ $
  float(STRMID(endTime,14,2))*60.+float(STRMID(endTime,17,2))
  startFile = where(time_sec ge startTimeSec)
  startFile = startFile[0]
  endFile = where(time_sec le endTimeSec)
  endFile = endFile[n_elements(endFile)-1]
  header_time = header_time[startFile:endFile]
  header_day = header_day[startFile:endFile]
  time_sec=time_sec[startFile:endFile]
  fileNameList=fileNameList[startFile:endFile]
  
  noTimeInput:
  
  ;------------------------ Read MVI data -------------------------------
  n_images = n_elements(time_sec)
  data_cube = fltarr(file_hdr.nx, file_hdr.ny, n_images)
  for i=startfile, endfile-1 do begin
    data_cube[*,*,i-startfile] = imgs[i]
  endfor
  CLOSE, /all ; Close 'lu' opened with OPENR
  
  ;-------------------- If Running Difference ---------------------------
  if ~exist(runDiff) then runDiff=3
  data_cube = data_cube-shift(data_cube,0,0,runDiff)
  image_type = 'RUNNING DIFF'
  
  ;-------------------- Reset some variables ----------------------------
  Restart0: ; After click in Reset (Distance x Time Plot)
  if exist(previousTime) then a0 = temporary(previousTime)
  if exist(angwidthset) then a0 = temporary(angwidthset)
  
  ReStart:  ; After click in change distance from the Sun (Angle x Time Plot)
  
  ;---------------- Check some input parameters -------------------------
  if ~keyword_set(angRange) then  angRange = [0,360]
  if ~keyword_set(xcen) then xcen = file_hdr.SUNXCEN
  if ~keyword_set(ycen) then ycen = file_hdr.SUNYCEN
  if ~keyword_set(sampleTime) then stepSize = 900. else stepSize = sampleTime
  
  ;----------------- Determinate Field of View in pixels ------------
  fovInPixel = (file_hdr.NX - xCen) < xCen < (file_hdr.NY - yCen) < ycen
  
  ;----------------- Circumference equation -----------------------------
  if ~keyword_set(cdelt1) then cdelt1 = file_hdr.SEC_PIX   ; Plate Scale
  
  if detector eq 'C2' then begin
  yymmdd=strmid(header_day[0],2,2)+strmid(header_day[0],5,2)+strmid(header_day[0],8,2)
  solar_ephem, yymmdd, radius=rSun_solar_ephem, /soho
  rSun=rSun_solar_ephem*60.*60.
  endif else begin
  rsun = header[0].RSUN      ; in arcsec (for COR2)
  endelse
  
  rSunInPx = rsun/cdelt1
  radius = inRadius*rSunInPx
  
  ; Check the limit of the radius
  if inradius lt 1 then begin
    print,'Min. radius input is: ', STRTRIM(STRING(1),1)
    inRadius = 1
    radius = inRadius*rSunInPx
  endif
  
  if radius gt fovInPixel then begin
    print,'Max. radius input is: ', STRTRIM(STRING(fovInPixel/rSunInPx,FORMAT='(f5.2)'),1)
    inRadius = inRadius-.5
    radius = inRadius*rSunInPx
  endif
  
  theta = (findgen(angRange[1]-angRange[0])+angRange[0])*!dtor+90.*!dtor
  x = radius*cos(theta)+xcen
  y = radius*sin(theta)+ycen
  
  ;------------------- Create the Jmap ----------------------------------
  n_images = (size(data_cube,/dim))[2]
  jmap = fltarr(n_images,n_elements(x))
  vecjmap = fltarr(n_elements(x))
  for j = 0,n_images-1 do begin
    for i = 0,n_elements(x)-1 do begin
      vecJmap[i] = float(data_cube[round(x[i]),round(y[i]),j])
    endfor
    jmap[j,*] = vecJmap
  endfor
  
  ;----------------- Create Time axis -----------------------------------
  totaltime = (time_sec[n_elements(time_sec)-1]-time_sec[0])/60.            ; in minutes
  xaxis = float([0,totaltime/5.,totaltime/5.*2,totaltime/5.*3,totaltime/5.*4,totaltime]) ; X-axis for HTmap plot
  
  ;----------------- Plot COR2 image------------------------------------
  window, 0, xs=wsize_x, ys=wsize_y
  plot_image,hist_equal(data_cube[*,*,8],per=1),back=255,color=0
  
  ; Plot Circumference and Angles 
  oplot,x,y,color=255
  xyouts,x[0]+5,y[0],STRTRIM(STRING(angRange[0]),1),color=255
  if angRange[1]-angRange[0] ge 360 then begin
    xyouts,x[n_elements(x)/4.]-5,y[n_elements(y)/4.],STRTRIM(STRING(angRange[1]/4),1),color=255, ALIGNMENT=1
  endif else begin
    xyouts,x[n_elements(x)-1.]-5.,y[n_elements(y)-1.],STRTRIM(STRING(angRange[1]),1),color=255, ALIGNMENT=1
  endelse
  
  ;--------------------- Fix X-axis w/ gaps -----------------------
  npixel = 1
  time_sec_old = time_sec[*,*]-time_sec[0,0]
  time_sec[*,*] = round((time_sec[*,*]-time_sec[0,0])/stepsize)*stepsize
  RMSD = sqrt(total((time_sec_old-time_sec)^2.)/n_images)
  
  xaxis_time = time_sec[n_images-1]         ; Time in seconds
  xaxis_size = xaxis_time/stepsize+1          ; Quantity of columns in the new J-map
  xaxis_position = time_sec/stepsize        ; Position of each columm of the old J-map in the new J-map
  
  sizejmap = size(jmap)
  slit_length = sizejmap[2]
  
  jmap_new0 = make_array(xaxis_size*npixel,slit_length)
  jmap_new0 = make_array(xaxis_size*npixel,360)
  
  for i = 0,n_images-1 do begin
    jmap_new0[xaxis_position[i]*npixel:xaxis_position[i]*npixel+npixel-1,angRange[0]:angRange[1]-1] $
    = jmap[i*npixel:i*npixel+npixel-1,*]
  endfor
  
  
  
  
  ; If you select certain inRadius close to Sun it will avoid the error in hist_equal when min=max
  if max(jmap_new0) eq min(jmap_new0) then begin
    dum_jmap_new0=jmap_new0
  endif else begin
    dum_jmap_new0 = hist_equal(jmap_new0,per=1)
  endelse
  
  ; If the Angle range is not 360, it will clean the area outside the main region
  if angRange[1]-angRange[0] ne 360 then begin
    dum_jmap_new0[*,0:angRange[0]]=255
    dum_jmap_new0[*,angRange[1]-1:359]=255
  endif
  
  window, 1, xs=wsize_x, ys=wsize_y
  itab0=41   ; Usually the last color index is 40, so it will change the Color Table index to 41.
  fg = 1b
  
  ; It will work with the Angle x Time and COR2 frames (Window 1 and 2)
  while fg do begin
  
    ; Load the itab0 only if it was manually defined
    ; If flagContrast0 is ON, means you want to call the XCOLORS (click in the right side of the plot)
    if exist(flagContrast0) then begin
      if ~flagContrast0 then loadct, itab0,/silent, file = '/Volumes/disk1/default.tbl'  
    endif
    
    wset, 1
    plot_image, dum_jmap_new0, ytitle='Angle (degree)', $
    xtitle='Elapsed time (min) since ' + header_day[0] + ' ' + strmid(header_time[0],0,8), $
      xticks=5, xtickname=STRTRIM(string(long(xaxis)),1),/smooth,back=255,color=0,scale=[3,1], $
      title = 'Distance from the Sun (Solar radii): ' + STRTRIM(string(inRadius,format='(f6.2)'),1), $
      yrange=[angrange[0],angrange[1]]
    loadct,0,/silent
    
    ; If you click in the right side of the plot (Open Color Table box)
    if keyword_set(flagContrast0) then begin
      flagContrast0=0b
      pcolor=!p.color
      !p.color=0
      XColors, ColorInfo=colorInfoPtr, NotifyPro='Refresh_Image', Image=dum_jmap_new0, WID=!D.Window, /block, $
        xtitle='Elapsed time (min) since ' + header_day[0] + ' ' + strmid(header_time[0],0,8), $
        ytitle='Angle (degree)', $
        xticks=5, xtickname=STRTRIM(string(long(xaxis)),1),/smooth,back=255,scale=[3,1], $
        imgtitle = 'Distance from the Sun (Solar radii): ' + STRTRIM(string(inRadius,format='(f6.2)'),1), $
        yrange=[angrange[0],angrange[1]]
      !p.color=pcolor
      ; The 'title=' is not in the '_Extra' because 'title' is also an input parameter for the Color Table box.
      MODIFYCT, itab0,'Xcolors set '+ STRTRIM(string(itab0),1), colorInfoPtr.R, colorInfoPtr.G, colorInfoPtr.B, $
                          file = '/Volumes/disk1/default.tbl' ; Save the new Color Table in itab0
    endif
    itab1=itab0+1
    loadct,0,/silent
    
    if exist(previousTime) then begin
      wset,1 & plots,previousTime,previousPos,psym=6,color=0,thick=2,symsize=.5
    endif
    if exist(angWidthSet) then begin
      if n_elements(angWidthSet) eq 1 then begin
        plots,timeAngWidthSet,angwidthSet,psym=1,color=0,thick=2
      endif else begin
        plots,timeAngWidthSet[0],angwidthSet[0],psym=1,color=0,thick=2
        plots,timeAngWidthSet[1],angwidthSet[1],psym=1,color=0,thick=2
      endelse
    endif
    
    ;----------------------- Get angle of the CME ---------------------
    wset,1
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; wait,.1
    cursor, curTime, curAngle
    if !mouse.button eq 2 then cursor, curTime, curAngle, /up
    ;cursor, curTime, curAngle,/down
    
    ; Increase and decrease the inRadius (distance from the Sun)
    ; when you click in the Left side of the Angle x Time Plot
    if curTime lt -10 then begin
      flagReStart=1
      if !mouse.button eq 1 then begin
        inRadius=inRadius+.5
        goto,ReStart
      endif
      if !mouse.button eq 4 then begin
        inRadius=inRadius-.5
        goto,ReStart
      endif
    endif
    
    ; When you click in the Right side of the Angle x Time plot,
    ; It set the flagContrast0 to open the Color Table box
    if curTime ge (size(dum_jmap_new0,/dim))[0]*3-3 then begin
      curTime = (size(dum_jmap_new0,/dim))[0]*3-3 ; Just to avoid error
      flagContrast0=1b
    endif
    
    ; If Left Click on the Angle x Time Plot
    if !mouse.button eq 1 then begin
      plots,curTime,curAngle,psym=1,color=255
      
      print,'Angle: ' + STRTRIM(string(long(curAngle)),1)+'º'
      
      ;-------------------- Create 1st order function -------------------
      theta = (curAngle+90.)*!dtor
      radius = findgen(fovInPixel)
      x1 = radius*cos(theta)+xcen
      y1 = radius*sin(theta)+ycen
      
      index = where(round(curtime/3.) le xaxis_position )
      wset,0
      ; Correspondent image to the time determined in the Angle/Time Map
      
      if ~keyword_set(flagContrast0) then begin
        plot_image,hist_equal(data_cube[*,*,index[0]],per=.1),back=255,color=0
        
        oplot,x,y,color=255
        xyouts,x[0]+5,y[0],STRTRIM(STRING(angRange[0]),1),color=255
        if angRange[1]-angRange[0] ge 360 then begin
          xyouts,x[n_elements(x)/4.]-5,y[n_elements(y)/4.],STRTRIM(STRING(angRange[1]/4),1),color=255, ALIGNMENT=1
        endif else begin
          xyouts,x[n_elements(x)-1.]-5.,y[n_elements(y)-1.],STRTRIM(STRING(angRange[1]),1),color=255, ALIGNMENT=1
        endelse
        oplot,x1,y1,color=255
      endif
      
      ; Plot the angular width
      if exist(angwidthSet) then begin
        if size([angwidthSet],/dim) eq 1 then begin
          aw_theta = (angwidthSet+90.)*!dtor
          x_aw = radius*cos(aw_theta)+xcen
          y_aw = radius*sin(aw_theta)+ycen
          oplot,x_aw,y_aw,color=180
        endif
        if size([angwidthSet],/dim) eq 2 then begin
          aw_theta = (angwidthSet[0]+90.)*!dtor
          x_aw = radius*cos(aw_theta)+xcen
          y_aw = radius*sin(aw_theta)+ycen
          oplot,x_aw,y_aw,color=180
          aw_theta=(angwidthSet[1]+90.)*!dtor
          x_aw = radius*cos(aw_theta)+xcen
          y_aw = radius*sin(aw_theta)+ycen
          oplot,x_aw,y_aw,color=180
        endif
      endif
      curtime_previous = CurTime
      curangle_previous = CurAngle
    endif
    
    ; If right click then finish the loop
    if !mouse.button EQ 4 THEN fg = 0b
    
    ; If middle click select Pre-event and Angles (Width)
    if !mouse.button eq 2 then begin
      if ~exist(previousTime) then begin
        previousTime = curTime
        previousPos = curAngle
      endif else begin
        angwidth = CurAngle
        angwidthSet = exist(angwidthSet) ? [angwidthSet,angwidth] : angwidth
        timeAngWidth = curTime
        timeAngWidthSet = exist(timeAngWidthSet) ? [timeAngWidthSet,curTime] : curTime
      endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      wait,.1     ; Necessary if use ctrl+click instead of middle button
    endif
    
    !mouse.button = 0 ; Reset click
    
    ; It fix the curTime and CurAngle for you to choose the Position Angle
    if keyword_set(curtime_previous) then begin
      curTime = curTime_previous
      curAngle = curAngle_previous
    endif
  endwhile
  
  ;------------ Determine total distance and Create Y-axis ----------
  totalDistInPixel = sqrt((x1[fovInPixel-1.]-x1[0])^2.+(y1[fovInPixel-1.]-y1[0])^2.)
  totalDistInArcsec = totalDistInPixel*cdelt1
  
  rSunInKm = onersun()
  DistSunEarth = rSunInKm/tan(rsun/3600.*!dtor) ; in km
  totalDistInKm = DistSunEarth*tan(totalDistInArcsec/3600.*!dtor)
  totalDistInGm = totalDistInKm/(10.^6.)
  
  totalDistInSolarRadii = totalDistInKm/rSunInKm
  yaxisSolarRadii = [0.,totalDistInSolarRadii/5.,totalDistInSolarRadii/5.*2.,totalDistInSolarRadii/5.*3.,totalDistInSolarRadii/5.*4.,totalDistInSolarRadii]
  
  yaxis = float([0.,totalDistInGm/5.,totalDistInGm/5.*2.,totalDistInGm/5.*3.,totalDistInGm/5.*4.,totalDistInGm])
  
  ;--------- Make the Jmap in the direction of the CME propagation -----------------
  jmap = fltarr(n_images,n_elements(x1))
  vecjmap = fltarr(n_elements(x1))
  for j = 0,n_images-1 do begin
    for i = 0,n_elements(x1)-1 do begin
      vecJmap[i] = float(data_cube[round(x1[i]),round(y1[i]),j])
    endfor
    jmap[j,*] = vecJmap
  endfor
  
  ;--------------------- Include gaps in the Jmaps -----------------------
  npixel = 1
  time_sec_old = time_sec[*,*]-time_sec[0,0]
  time_sec[*,*] = round((time_sec[*,*]-time_sec[0,0])/stepsize)*stepsize
  RMSD = sqrt(total((time_sec_old-time_sec)^2.)/n_images)
  
  xaxis_time = time_sec[n_images-1]         ; Time in seconds
  xaxis_size = xaxis_time/stepsize+1          ; Quantity of columns in the new J-map
  xaxis_position = time_sec/stepsize        ; Position of each columm of the old J-map in the new J-map
  
  sizejmap = size(jmap)
  slit_length = sizejmap[2]
  
  jmap_new = make_array(xaxis_size*npixel,slit_length)
  
  jmapinput = jmap ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Modified 12/18
  for i = 0,n_images-1 do begin
    jmap_new[xaxis_position[i]*npixel:xaxis_position[i]*npixel+npixel-1,*] = jmapinput[i*npixel:i*npixel+npixel-1,*]
  endfor
  
  ;--------------------------- Plot COR2 image ------------------------------
  index = where(round(curtime/3) le xaxis_position )
  
  window, 0, xs  = wsize_x, ys = wsize_y   ; este plot quedo redundante
  plot_image,hist_equal(data_cube[*,*,index[0]],per=.1),back=255,color=0
  
  oplot,x,y,color = 255
  xyouts,x[0]+5,y[0],STRTRIM(STRING(angRange[0]),1),color=255
  
  if angRange[1]-angRange[0] ge 360 then begin
    xyouts,x[n_elements(x)/4.]-5,y[n_elements(y)/4.],STRTRIM(STRING(angRange[1]/4),1),color=255, ALIGNMENT=1
  endif else begin
    xyouts,x[n_elements(x)-1.]-5.,y[n_elements(y)-1.],STRTRIM(STRING(angRange[1]),1),color=255, ALIGNMENT=1
  endelse
  oplot,x1,y1,color=255
  
  ; Plot the angular width
  if exist(angwidthSet) then begin
    if size([angwidthSet],/dim) eq 1 then begin
      aw_theta = (angwidthSet+90.)*!dtor   ; angular_width_theta
      x_aw = radius*cos(aw_theta)+xcen
      y_aw = radius*sin(aw_theta)+ycen
      oplot,x_aw,y_aw,color=180
    endif
    if size([angwidthSet],/dim) eq 2 then begin
      aw_theta = (angwidthSet[0]+90.)*!dtor
      x_aw = radius*cos(aw_theta)+xcen
      y_aw = radius*sin(aw_theta)+ycen
      oplot,x_aw,y_aw,color=180
      aw_theta = (angwidthSet[1]+90.)*!dtor
      x_aw = radius*cos(aw_theta)+xcen
      y_aw = radius*sin(aw_theta)+ycen
      oplot,x_aw,y_aw,color = 180
    endif
  endif
  
  ; If you want compare the time of the two HTmap.
  if keyword_set(compare) then begin
    ref_size=size(jmap_new,/dim)
    comp_jmap_new=congrid(hist_equal(jmap_new,per=1), 3*ref_size[0],ref_size[1], /cubic)
    comp_dum_jmap_new0=congrid(dum_jmap_new0, 3*ref_size[0],ref_size[1], /cubic)
    comp_img=bytarr(3*ref_size[0],2*ref_size[1])
    comp_img[*,0:ref_size[1]-1]=comp_dum_jmap_new0
    comp_img[*,ref_size[1]:2*ref_size[1]-1]=comp_jmap_new
    comp_img[*,ref_size[1]]=bytarr(3*ref_size[0])
    window, 8, xs=450, ys=650
    plot_image,comp_img,color=0,back=255
  endif
  

 ;----------------  While -------------- 
 window, 2, xs =wsize_x, ys = wsize_y
  window, 4;, xs =wsize_x, ys = wsize_y
  fg = 1b
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;loadct,0,/silent
  while fg do begin
  
    if exist(flagContrast1) then begin
      if ~flagContrast1 then loadct, itab1,/silent, file= '/Volumes/Disk1/default1.tbl' 
    endif
    wset, 4
    plot_image,hist_equal(jmap_new,per=1),ytitle='Distance from the center of the sun (1M km)', $
      xtitle='Elapsed time (min) since ' + header_day[0] + ' ' + strmid(header_time[0],0,8), $
      xticks=5, xtickname=STRTRIM(string(long(xaxis)),1), $
      yticks=5, ytickname=STRTRIM(string(float(yaxis),FORMAT='(f5.2)'),1), $
      /smooth,back=255,color=0,scale=[5,2], $
      xmargin=[3.3,3.3],ystyle=1+8
    axis,yaxis=1,yticks=5,ytickn=STRTRIM(string(float(yaxisSolarRadii),FORMAT='(f5.2)'),1), $
      ystyle=1,color=0,ytitle= 'Distance (Solar Radius)'
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; loadct,0,/silent
    
    if keyword_set(flagContrast1) then begin
      flagContrast1=0b
      pcolor=!p.color
      !p.color=0
      XColors, ColorInfo=colorInfoPtr, NotifyPro='Refresh_Image', Image=hist_equal(jmap_new,per=1), WID=!D.Window, /block, $
      xtitle='Elapsed time (min) since ' + header_day[0] + ' ' + strmid(header_time[0],0,8), $
      ytitle='Distance from the center of the sun (1M km)', $
        xticks=5, xtickname=STRTRIM(string(long(xaxis)),1), $
        yticks=5, ytickname=STRTRIM(string(float(yaxis),FORMAT='(f5.2)'),1), $
        /smooth,back=255,scale=[5,2], $
        xmargin=[3.3,3.3],ystyle=1+8,yaxisSolarRadii=yaxisSolarRadii
      !p.color=pcolor
      MODIFYCT, itab1,'Xcolors set '+ STRTRIM(string(itab0),1), colorInfoPtr.R, colorInfoPtr.G, colorInfoPtr.B, $
                           file = '/Volumes/Disk1/default1.tbl'       
    endif
    
    if exist(xx0) then begin
      for j = 0,(size([A],/dim))[0]-1 do begin
        xyouts,xx0MaxSet[j]*5,yy0MaxSet[j]*2,STRTRIM(string(long(j+1)),1)
      endfor
    endif
    
    wait,.09
    cursor, curTime1, curDist
    
    if (CurTime1 lt -40) or (curDist lt 0) then goto, ReStart0
    
    ; If click outside the image (between 0 and -40). It is a margin for error.
    if curTime1 lt 0 then curTime1 = 0
    
    ; If click outside the image (right side) fix the curTime1
    if curTime1 ge (size(jmap_new,/dim))[0]*5-5 then begin
      curTime1 = (size(jmap_new,/dim))[0]*5-5
      flagContrast1=1b
    endif
    
    ; Now we need the position in jmapfilter, that are equal the number of the image
    curDist = curDist/2
    
    wset, 2   ; It is redundant
    
    again:
    idx = where(xaxis_position eq  round(curtime1/5.), ct)   ; ct for use in case it does not find.
    if ct eq 0 then begin
      curtime1 = curtime1+5
      goto, again
    endif
    
    ;----------------- Circumference equation -----------------------------
    theta = ((findgen(angRange[1]-angRange[0])+angRange[0])+90.)*!dtor
    radio = curDist
    x = radio*cos(theta)+xcen
    y = radio*sin(theta)+ycen
    
    sz = size(data_cube,/dim)
    img = bytarr(sz[0],sz[1],3)
    ;img[*,*,0] = hist_equal(data_cube[*,*,idx[0]-1], per=.1)
    img[*,*,1] = hist_equal(data_cube[*,*,idx[0]], per=.1)
    ;img[*,*,2] = hist_equal(data_cube[*,*,idx[0]+1], per=.1)
    if ~keyword_set(flagContrast1) then begin
      plot_image, hist_equal(img[*,*,1], per=.1),back=255,color=0
      
      oplot,x,y,color=255
      xyouts,x[0]+5,y[0],STRTRIM(STRING(angRange[0]),1),color=255
      if angRange[1]-angRange[0] ge 360 then begin
        xyouts,x[n_elements(x)/4.]-5,y[n_elements(y)/4.],STRTRIM(STRING(angRange[1]/4),1),color=255, ALIGNMENT=1
      endif else begin
        xyouts,x[n_elements(x)-1.]-5.,y[n_elements(y)-1.],STRTRIM(STRING(angRange[1]),1),color=255, ALIGNMENT=1
      endelse
      oplot,x1,y1,color=255
    endif
    ; If middle click
    if !mouse.button eq 2 then begin
      if exist(flagContrast1) then loadct,itab1,/silent, file= '/Volumes/Disk1/default1.tbl' 
      wset, 4   ; If I don't re-plot the image, it screws things up.
      plot_image,hist_equal(jmap_new,per=1),xtitle='Elapsed time (min) since ' + header_day[0] + ' ' + strmid(header_time[0],0,8), $
      ytitle='Distance from the center of the sun (1M km)', $
        xticks=5, xtickname=STRTRIM(string(long(xaxis)),1), $
        yticks=5, ytickname=STRTRIM(string(float(yaxis),FORMAT='(f5.2)'),1),$
        /smooth,back=255,color=0,scale=[5,2], $
        xmargin=[3.3,3.3],ystyle=1+8
      axis,yaxis=1,yticks=5,ytickn=STRTRIM(string(float(yaxisSolarRadii),FORMAT='(f5.2)'),1), $
        ystyle=1,color=0,ytitle= 'Distance (Solar Radius)'
      loadct,0,/silent
      clickjmap_cor2_test,jmap_new,spfile=spfile
      remark='' & read, remark, prompt='REMARK: '
      remarkSet = exist(remarkSet) ? [remarkset,remark] : remark
      
      if (size(spfile))[0] ne 1 then newspeed=fltarr((size(spfile,/dim))[0]+4,(size(spfile,/dim))[1]) else $
        newspeed = fltarr((size(spfile,/dim))[0]+4)
           
      x = spfile[0:1,0]
      for k = 1,(size(spfile))[2]-1 do x=[x,spfile[0:1,k]]
      y = spfile[2:3,0]
      for k = 1,(size(spfile))[2]-1 do y=[y,spfile[2:3,k]]
      
      coefs = poly_fit(x,y,2,yband=yband, sigma=sigma)
      xx = findgen(1000)*(max(x)-min(x))/999.+min(x)
      yy = coefs[2]*xx^2.+coefs[1]*xx+coefs[0]
      
      xx0 = xx/5.
      srtotal = totalDistInArcsec 
      yy0 = yy/2.
      
      ;--- To be used in the Set of COR2 images
      Seq_timeInMin = exist(Seq_timeInMin) ? [[Seq_timeInMin],[(xx0-1.)*stepSize/60.]] : (xx0-1.)*stepSize/60.  ; In minutes
      Seq_radius = exist(Seq_radius) ? [[Seq_radius],[yy0]] : yy0           ; In pixels
      ;---
      
      jmap_new1 = jmap_new
      jmap_new1[round(xx0),round(yy0)] = min(jmap_new1)
      
      ;------------ Poly_fit real data ------------------
      xreal = round(spfile[0:1,0]/5.)
      for k = 1,(size(spfile))[2]-1 do xreal=[xreal,round(spfile[0:1,k]/5)]
      yreal = round(spfile[2:3,0]/2.)
      for k = 1,(size(spfile))[2]-1 do yreal=[yreal,round(spfile[2:3,k]/2.)]
      
      xreal = xreal*totaltime*60/(size(jmap_new1,/dim))[0]
      yreal = yreal*totalDistInKm/(size(jmap_new1,/dim))[1]
      
      coefs = poly_fit(xreal,yreal,2,yband=yband, sigma=sigma)
      xxreal = findgen(1000)*(max(xreal)-min(xreal))/999.+min(xreal)
      yyreal = coefs[2]*xxreal^2.+coefs[1]*xxreal+coefs[0]
      
      coef1stOrder = poly_fit(xreal,yreal,1,yband=yband1storder, sigma=sigma1stOrder)
      
      setSpeed1stOrder = exist(Setspeed1stOrder) ? [Setspeed1stOrder,coef1stOrder[1]] : coef1stOrder[1]
      setSigma1stOrder = exist(setSigma1stOrder) ? [setSigma1stOrder,sigma1stOrder[1]] : sigma1stOrder[1]
      
      ; Coefficients of the 2nd Order equation
      A = exist(A) ? [A,coefs[2]] : coefs[2]
      B = exist(B) ? [B,coefs[1]] : coefs[1]
      C = exist(C) ? [C,coefs[0]] : coefs[0]
      Ea = exist(Ea) ? [Ea,sigma[2]/sqrt(n_elements(xreal))] : sigma[2]/sqrt(n_elements(xreal))
      Eb = exist(Eb) ? [Eb,sigma[1]/sqrt(n_elements(xreal))] : sigma[1]/sqrt(n_elements(xreal))
      Ec = exist(Ec) ? [Ec,sigma[0]/sqrt(n_elements(xreal))] : sigma[0]/sqrt(n_elements(xreal))
      xStart = exist(xStart) ? [xStart,xxreal[0]] : xxreal[0]
      xEnd = exist(xEnd) ? [xEnd,xxreal[(size(xxreal,/dim))[0]-1]] : xxreal[(size(xxreal,/dim))[0]-1]
      yStart = exist(yStart) ? [yStart,yyreal[0]] : yyreal[0]
      yEnd = exist(yEnd) ? [yEnd,yyreal[(size(yyreal,/dim))[0]-1]] : yyreal[(size(yyreal,/dim))[0]-1]
      
      ; Create the set of clicks.
      ; It's so extensive because you can create sets with different sizes.
      ; It means, different vectors, that we still don't know the size, in the same array.
      if exist(xAxisClickSet) then begin
        n_xreal = n_elements(xreal)
        n_xAxisClickSet = (size(xAxisClickSet))[0] eq 1 ?  [size(xAxisClickSet,/dim),1] : size(xAxisClickSet,/dim)
        increase = abs(n_xAxisClickSet[0]-n_xreal)
        SetNpoints = [SetNpoints,n_xreal]
        if n_xAxisClickSet[0] gt n_xreal then begin
          xAxisClickSet = [[xAxisClickSet],[reform(xreal),fltarr(increase)]]
          yAxisClickSet = [[yAxisClickSet],[reform(yreal),fltarr(increase)]]
          xclickpx = [[xclickpx],[reform(round(x/5.)),fltarr(increase)]]
          yclickpx = [[yclickpx],[reform(round(y/2.)),fltarr(increase)]]
        endif else begin
          if n_xAxisClickSet[0] lt n_xreal then begin
            xAxisClickSet = [[xAxisClickSet,fltarr(increase,n_xAxisClickSet[1])],[reform(xreal)]]
            yAxisClickSet = [[yAxisClickSet,fltarr(increase,n_xAxisClickSet[1])],[reform(yreal)]]
            xclickpx = [[xclickpx,fltarr(increase,n_xAxisClickSet[1])],[reform(round(x/5.))]]
            yclickpx = [[yclickpx,fltarr(increase,n_xAxisClickSet[1])],[reform(round(y/2.))]]
          endif else begin
            xAxisClickSet = [[xAxisClickSet],[reform(xreal)]]
            yAxisClickSet = [[yAxisClickSet],[reform(yreal)]]
            xclickpx = [[xclickpx],[reform(round(x/5.))]]
            yclickpx = [[yclickpx],[reform(round(y/2.))]]
          endelse
        endelse
      endif else begin
        xAxisClickSet = reform(xreal)
        yAxisClickSet = reform(yreal)
        xclickpx = reform(round(x/5.))                  ; Clicked points in Pixels
        yclickpx = reform(round(y/2.))
        SetNpoints = n_elements(xAxisClickSet)
      endelse
      
      if exist(xx0) then begin
        xx0MaxSet = exist(xx0MaxSet) ? [xx0MaxSet,xx0[(size(xx0,/dim))[0]-1]] : xx0[(size(xx0,/dim))[0]-1]
        yy0MaxSet = exist(yy0MaxSet) ? [yy0MaxSet,yy0[(size(yy0,/dim))[0]-1]] : yy0[(size(yy0,/dim))[0]-1]
      endif
      speedEquation2nd = STRTRIM(string(2*A,FORMAT='(f9.3)'),1) + ' +- ' + STRTRIM(string(2*Ea,FORMAT='(f9.3)'),1)  + ')*t(s) + (' + STRTRIM(string(B,FORMAT='(f11.3)'),1) + ' +- ' + STRTRIM(string(Eb,FORMAT='(f9.3)'),1) + ')'
      print,'Speed 2nd: ' + 'V(km/s) = (' + speedEquation2nd
      !mouse.button = 0
      jmap_new = jmap_new1
    endif
    if !mouse.button EQ 4 THEN fg = 0b
    !mouse.button=0
  endwhile
  
  ;---------------------- Speed Plot ---------------------
  if ~exist(seq_timeInMin) then goto,noSpeed
  
  loadct,5,/silent
  totalDistanceSpeedGraphic = 10.*10.^6./rSuninkm
  offsetDistance = 2.*10.^6./rSuninkm
  
  ;--------------- X-axis for the Speed Plot ---------
  xaxisSolarRadii = [0,3,6,9,12,15] ; These numbers are the sticks for the X-axis
  if detector eq 'C2' then xaxisSolarRadii = [2,3,4,5,6,7]
  xaxisSolarRadii2 = round(floor(xaxisSolarRadii))
  SpeedGraph_xtickv = (xaxisSolarRadii2*rSunInKm/(10.^6.))
  
  ;---------------- Determine the yrange ----------------------------
  for i = 0,(size([A],/dim))[0]-1 do begin
    timeAxis = findgen(1000)/999.*(xEnd[i]-xStart[i])+xStart[i]
    speed_2nd = 2*A[i]*timeAxis+B[i]
    max_speed_2nd = exist(max_speed_2nd) ? max(speed_2nd) > max_speed_2nd : max(speed_2nd)
  endfor
  
  ;-------------------- Speed x Distance Plot -------------------------
  for i = 0,(size([A],/dim))[0]-1 do begin
    timeAxis = findgen(1000)/999.*(xEnd[i]-xStart[i])+xStart[i]
    DistAxis = A[i]*timeAxis^2+B[i]*timeAxis+C[i]
    speed_2nd = 2*A[i]*timeAxis+B[i]
    errSup = 2*(A[i]+Ea[i])*timeAxis+(B[i]+Eb[i])
    errinf = 2*(A[i]-Ea[i])*timeAxis+(B[i]-Eb[i])
    if i eq 0 then begin
      window, 5
      plot, distAxis/10.^6., speed_2nd, thick=2, $
        xtitle = 'Distance from the center of the sun (1M km)', ytitle='Speed (km/s)',$
        background=255,color=0, $
        xstyle=8,ymargin=[3.3,3.3],xrange=detector eq 'C2'?[1,5]:[0,12],yrange=[0,max_speed_2nd]
      axis,xaxis=1,xticks=detector eq 'C2'?5:5,xtickn=STRTRIM(string(xaxisSolarRadii2),1), $
        xstyle=1,color=0,xtitle= 'Distance (Solar Radius)',xtickv=SpeedGraph_xtickv,xminor=3
      oplot, distAxis/10.^6., errSup, linestyle=2,color=0
      oplot, distAxis/10.^6., errinf, linestyle=2,color=0
      xyouts,.95,.14, /normal, alignment=1.,color=0, $
        'Speed 2nd: V(km/s) = (' + STRTRIM(string(2*A[i],FORMAT='(f9.4)'),1) + ' +- ' + STRTRIM(string(2*Ea[i],FORMAT='(f9.4)'),1)  + ')*t(s) + (' + STRTRIM(string(B[i],FORMAT='(f11.3)'),1) + ' +- ' + STRTRIM(string(Eb[i],FORMAT='(f9.3)'),1) + ')'
      xyouts,distAxis[n_elements(distAxis)-1]/10.^6.*1.005,speed_2nd[n_elements(errsup)-1]*0.987,'1',color=0
    endif else begin
      oplot, distAxis/10.^6., speed_2nd, thick=2,color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
      oplot, distAxis/10.^6., errSup, linestyle=2,color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
      oplot, distAxis/10.^6., errinf, linestyle=2,color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
      xyouts,distAxis[n_elements(distAxis)-1]/10.^6.*1.005,speed_2nd[n_elements(errsup)-1]*0.987,STRTRIM(string(i+1),1),color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
      xyouts,.95,.14+i*.05, /normal, alignment=1.,color=!D.TABLE_SIZE*i/(size(A,/dim))[0], $
        'Speed 2nd: V(km/s) = (' + STRTRIM(string(2*A[i],FORMAT='(f9.4)'),1) + ' +- ' + STRTRIM(string(2*Ea[i],FORMAT='(f9.4)'),1)  + ')*t(s) + (' + STRTRIM(string(B[i],FORMAT='(f11.3)'),1) + ' +- ' + STRTRIM(string(Eb[i],FORMAT='(f9.3)'),1) + ')'
    endelse
  endfor
  xyouts,.95,.14+i*.05 , /normal, alignment=1.,color=0, 'Position angle (deg) = ' + STRTRIM(string(curAngle,FORMAT='(f7.2)'),1)
  
  ;--------------------- Distance x Time Plot -----------------
  for i = 0,(size([A],/dim))[0]-1 do begin
    timeAxis = findgen(1000)/999.*(xEnd[i]-xStart[i])+xStart[i]
    distAxis = A[i]*timeAxis^2+B[i]*timeAxis+C[i]
    if i eq 0 then begin
      window, 6
      plot,timeAxis/60.,distAxis/10.^6., thick=2, $
        xtitle='Elapsed time (min) since ' + header_day[0] + ' ' + strmid(header_time[0],0,8), $
        ytitle='Distance from the center of the sun (1M km)', $
        background=255,color=0, $
        ystyle=8,xmargin=[6,6],yrange=detector eq 'C2'?[1,5]:[0,12]
      axis,yaxis=1,yticks=detector eq 'C2'?5:5,ytickn=STRTRIM(string(xaxisSolarRadii2),1), $
        ystyle=1,color=0,ytitle= 'Distance (Solar Radius)',ytickv=SpeedGraph_xtickv,yminor=3
      plots,xAxisClickSet[*,i]/60.,yAxisClickSet[*,i]/10.^6.,psym=1,color=0
      
      xyouts,.92,.14, /normal, alignment=1.,color=0, $
        'D(km) = (' + STRTRIM(string(A[i],FORMAT='(f9.4)'),1) + ' +- ' + STRTRIM(string(Ea[i],FORMAT='(f9.4)'),1)  + ')*t(s)^2 + (' + $
        STRTRIM(string(B[i],FORMAT='(f11.3)'),1) + ' +- ' + STRTRIM(string(Eb[i],FORMAT='(f9.3)'),1) + ')*t(s) + (' + $
        STRTRIM(string(C[i],FORMAT='(f12.3)'),1) + ' +- ' + STRTRIM(string(Ec[i],FORMAT='(f12.3)'),1) + ')',charsize=.9
      xyouts,timeAxis[n_elements(timeAxis)-1]/60.*1.005,distAxis[n_elements(distAxis)-1]/10.^6.*0.987,'1',color=0
      
    endif else begin
      oplot, timeAxis/60., distAxis/10.^6., thick=2,color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
      plots,xAxisClickSet[*,i]/60.,yAxisClickSet[*,i]/10.^6.,psym=1,color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
      xyouts,.92,.14+i*.045, /normal, alignment=1.,color=!D.TABLE_SIZE*i/(size(A,/dim))[0], $
        'D(km) = (' + STRTRIM(string(A[i],FORMAT='(f9.4)'),1) + ' +- ' + STRTRIM(string(Ea[i],FORMAT='(f9.4)'),1)  + ')*t(s)^2 + (' + $
        STRTRIM(string(B[i],FORMAT='(f11.3)'),1) + ' +- ' + STRTRIM(string(Eb[i],FORMAT='(f9.3)'),1) + ')*t(s) + (' + $
        STRTRIM(string(B[i],FORMAT='(f11.3)'),1) + ' +- ' + STRTRIM(string(Eb[i],FORMAT='(f9.3)'),1) + ')',charsize=.9
      xyouts,timeAxis[n_elements(timeAxis)-1]/60.*1.005,distAxis[n_elements(distAxis)-1]/10.^6.*0.987,STRTRIM(string(i+1),1),color=!D.TABLE_SIZE*i/(size(A,/dim))[0]
    endelse
  endfor
  xyouts,.92,.14+i*.045 , /normal, alignment=1.,color=0, 'Position angle (deg) = ' + STRTRIM(string(curAngle,FORMAT='(f7.2)'),1),charsize=.9
  
  ;--------------------------- Set of COR2 images ---------------------
  
  ;------------ Time and image index --------------
  SeqImageInst = [.2,.5,.8] ;Percent of the total time to plot on the Set image.
  sz_seq_timeInMin = (size([seq_timeInMin]))[0] eq 1 ? [(size([seq_timeInMin],/dim))[0],1] :  size([seq_timeInMin],/dim)
  seq_timeSample = [seq_timeInMin[(sz_seq_timeInMin[0]-1)*SeqImageInst[0],0],seq_timeInMin[(sz_seq_timeInMin[0]-1)*SeqImageInst[1],0],seq_timeInMin[(sz_seq_timeInMin[0]-1)*SeqImageInst[2],0]]
  
  seq_imageIndex = [(where(time_sec ge seq_timesample[0]*60))[0],(where(time_sec ge seq_timesample[1]*60))[0],(where(time_sec ge seq_timesample[2]*60))[0]]
  
  if sz_seq_timeInMin[1] gt 1 then begin
    seq_timeIndex2 = fltarr(3,sz_seq_timeInMin[1])
    for j = 0,sz_seq_timeInMin[1]-1 do begin
      seq_timeIndex2[*,j] = [(where(time_sec[seq_imageIndex[0]] le seq_timeInMin[*,j]*60.))[0],(where(time_sec[seq_imageIndex[1]] le seq_timeInMin[*,j]*60.))[0],(where(time_sec[seq_imageIndex[2]] le seq_timeInMin[*,j]*60.))[0]]
    endfor
  endif
  
  ;------------ Using time index in the radius vector --------------
  sz_seq_radius = (size([seq_radius]))[0] eq 1 ? [(size([seq_radius],/dim))[0],1] :  size([seq_radius],/dim)
  
  if sz_seq_timeInMin[1] gt 1 then begin
    seq_radiusSample = fltarr(3,sz_seq_radius[1])
    for j = 0,sz_seq_radius[1]-1 do begin
      seq_radiusSample[*,j] = [seq_radius[seq_timeIndex2[0,j],j],seq_radius[seq_timeIndex2[1,j],j],seq_radius[seq_timeIndex2[2,j],j]]
    endfor
  endif else begin
    seq_radiusSample = [seq_radius[(sz_seq_timeInMin[0]-1)*SeqImageInst[0]],seq_radius[(sz_seq_timeInMin[0]-1)*SeqImageInst[1]],seq_radius[(sz_seq_timeInMin[0]-1)*SeqImageInst[2]]]
  endelse
  
  ;------------ Create the array with the set of the image --------------
  seq_imageSet = bytarr((size(data_cube,/dim))[0]*3,(size(data_cube,/dim))[1])
  seq_imageSet[0:(size(data_cube,/dim))[0]-1,*] = hist_equal(data_cube[*,*,seq_imageindex[0]]-data_cube[*,*,seq_imageindex[0]-3],per = 1)
  seq_imageSet[(size(data_cube,/dim))[0]:(size(data_cube,/dim))[0]*2-1,*] = hist_equal(data_cube[*,*,seq_imageindex[1]]-data_cube[*,*,seq_imageindex[1]-3],per=1)
  seq_imageSet[(size(data_cube,/dim))[0]*2:(size(data_cube,/dim))[0]*3-1,*] = hist_equal(data_cube[*,*,seq_imageindex[2]]-data_cube[*,*,seq_imageindex[2]-3],per=1)
  
  ;------------ Plot the set, degree and radius --------------
  loadct,0,/silent
  window,7,xs=1000,ys=400
  plot_image,seq_imageSet,back=255,color=0,xstyle=4+1,ystyle=4+1,xmargin=[1,1],ymargin=[1,1]
  curAngle = curangle+90; This way the 0deg start in the top (it's used only during the plot)
  
  circ_theta = (findgen(angRange[1]-angRange[0])+angRange[0])*!dtor+90.*!dtor
  circ_radio = inRadius*rSunInPx
  circ_x = circ_radio*cos(circ_theta)
  circ_y = circ_radio*sin(circ_theta)
  
  loadct,5,/silent ; Load the Color Table 5
  for j = 0,sz_seq_radius[1]-1 do begin
    for i = 0,2 do begin
      if j eq 0 then begin
        plots,radius*cos(curAngle*!dtor)+xcen+(size(data_cube,/dim))[0]*i,radius*sin(curAngle*!dtor)+ycen,color=255
        plots,seq_radiusSample[i,j]*cos(curAngle*!dtor)+xcen+(size(data_cube,/dim))[0]*i,seq_radiusSample[i,j]*sin(curAngle*!dtor)+ycen,psym=1,color=0,thick=2
        oplot,circ_x+xcen+(size(data_cube,/dim))[0]*i,circ_y+ycen,color=255
        xyouts,circ_x[0]+xcen+(size(data_cube,/dim))[0]*i+5,circ_y[0]+ycen,STRTRIM(STRING(angRange[0]),1),color=255
        xyouts,circ_x[n_elements(x)/4.]+xcen+(size(data_cube,/dim))[0]*i-5,circ_y[n_elements(y)/4.]+ycen,STRTRIM(STRING(angRange[1]/4),1),color=255, ALIGNMENT=1
      endif else begin
        plots,seq_radiusSample[i,j]*cos(curAngle*!dtor)+xcen+(size(data_cube,/dim))[0]*i,seq_radiusSample[i,j]*sin(curAngle*!dtor)+ycen,psym=1, $
          color=!D.TABLE_SIZE*j/sz_seq_radius[1],thick=2
        oplot,circ_x+xcen+(size(data_cube,/dim))[0]*i,circ_y+ycen,color=255
        xyouts,circ_x[0]+xcen+(size(data_cube,/dim))[0]*i+5,circ_y[0]+ycen,STRTRIM(STRING(angRange[0]),1),color=255
        xyouts,circ_x[n_elements(x)/4.]+xcen+(size(data_cube,/dim))[0]*i-5,circ_y[n_elements(y)/4.]+ycen,STRTRIM(STRING(angRange[1]/4),1),color=255, ALIGNMENT=1
      endelse
    endfor
  endfor
  
  curAngle = curangle-90 ; Reset curAngle after plot to the horizontal reference.
  loadct,0,/silent ; Reset Color Table
  
  ;------------------------ Save the images -----------------------
  wset,1
  fnameout1 ='HT_AngleTime' + '_' + STRMID(header_day[0],0,4) + STRMID(header_day[0],5,2) + STRMID(header_day[0],8,2) $
    + '_' + STRMID(header_time[0],0,2) + STRMID(header_time[0],3,2) + STRMID(header_time[0],6,2) + '.png'
  write_png, mvipath + fnameout1,tvrd()
  wset,4
  fnameout4 ='HT_DistanceTime' + '_' + STRMID(header_day[0],0,4) + STRMID(header_day[0],5,2) + STRMID(header_day[0],8,2) $
    + '_' + STRMID(header_time[0],0,2) + STRMID(header_time[0],3,2) + STRMID(header_time[0],6,2) + '.png'
  write_png, mvipath + fnameout4,tvrd()
  wset,5
  fnameout5 ='Plot_SpeedDistance' + '_' + STRMID(header_day[0],0,4) + STRMID(header_day[0],5,2) + STRMID(header_day[0],8,2) $
    + '_' + STRMID(header_time[0],0,2) + STRMID(header_time[0],3,2) + STRMID(header_time[0],6,2) + '.png'
  write_png, mvipath + fnameout5,tvrd(true=1)
  wset,6
  fnameout6 ='Plot_DistanceTime' + '_' + STRMID(header_day[0],0,4) + STRMID(header_day[0],5,2) + STRMID(header_day[0],8,2) $
    + '_' + STRMID(header_time[0],0,2) + STRMID(header_time[0],3,2) + STRMID(header_time[0],6,2) + '.png'
  write_png, mvipath + fnameout6,tvrd(true=1)
  wset,7
  fnameout7 ='ImageSequence' + '_' + STRMID(header_day[0],0,4) + STRMID(header_day[0],5,2) + STRMID(header_day[0],8,2) $
    + '_' + STRMID(header_time[0],0,2) + STRMID(header_time[0],3,2) + STRMID(header_time[0],6,2) + '.png'
  write_png, mvipath + fnameout7,tvrd(true=1)
  
  ;-------------- Record the catalog text ----------
  cor2record_test, header_day, header_time, time_sec0, xstart, xend, detector, image_type, curangle, $
    speedEquation2nd, xAxisClickSet, yAxisClickSet, xclickpx, yclickpx, setnpoints, cdelt1, rsun, $
    n_images, A, Ea, B, Eb, setSpeed1stOrder, setSigma1stOrder, fnameout1, fnameout4, fnameout5, $
    fnameout6, fnameout7,remarkSet, angwidthset, timeAngWidthSet, previousTime, stepSize, inRadius, $
    mviFileName, mvipath, fileNameList, xaxis_position
    
  noSpeed:
  stop
  
end