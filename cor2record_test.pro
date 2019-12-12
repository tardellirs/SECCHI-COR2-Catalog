pro cor2record_test, header_day, header_time, time_sec0, xstart, xend, detector, image_type, curangle, $
    speedEquation2nd, xAxisClickSet, yAxisClickSet, xclickpx, yclickpx, setnpoints, cdelt1, rsun, $
    n_images, A, Ea, B, Eb, setSpeed1stOrder, setSigma1stOrder, fnameout1, fnameout4, fnameout5, $
    fnameout6, fnameout7, remarkSet, angwidthset, timeAngWidthSet, previousTime, stepSize, inRadius, $
    mviFileName, mvipath, fileNameList, xaxis_position
    
  ;+
  ;
  ;   Name: cor2record_test
  ;
  ;   Purpose: Record informations in a Text File for the COR2 Catalog
  ;
  ;-
    
  path= mvipath
  filename='COR2CAT' + '_' + STRMID(header_day[0],0,4) + STRMID(header_day[0],5,2) + STRMID(header_day[0],8,2) $
    + '_' + STRMID(header_time[0],0,2) + STRMID(header_time[0],3,2) + STRMID(header_time[0],6,2) + '.txt'
  fid=1
  
  ;-- Determine PLOT_START_TIME:
  plot_start_time = header_day[0] + ' ' + STRMID(header_time[0],0,8)
  
  ;-- Determine PLOT_END_TIME:
  plot_end_time = header_day[n_elements(header_day)-1] + ' ' + STRMID(header_time[n_elements(header_time)-1],0,8)
  
  ;-- Determine DATE-OBS:
  time_sec=time_sec0[*,*]-time_sec0[0,0]
  date_obs=strarr(n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    idxDateObs=(where(xstart[i] le time_sec))[0]
    date_obs[i]=header_day[idxDateObs]
  endfor
  
  ;-- Determine TIME-OBS:
  time_obs=strarr(n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    time_obs_sec = xstart[i]+time_sec0[0]-long(STRMID(header_day[0],8,2))*60*60*24
    time_obs_hr = floor(time_obs_sec/(60.*60.))
    time_obs_min = floor(time_obs_sec/60.-time_obs_hr*60.)
    sTime_obs_sec = STRTRIM(string(time_obs_sec - time_obs_min*60. - time_obs_hr*60.*60,format='(i2.2)'),1)
    sTime_obs_min = STRTRIM(string(time_obs_min,format='(i2.2)'),1)
    if long(time_obs_hr) ge 24 then time_obs_hr = 0
    sTime_obs_hr = STRTRIM(string(time_obs_hr,format='(i2.2)'),1)
    time_obs[i] = sTime_obs_hr + ':' + sTime_obs_min + ':' + sTime_obs_sec
  endfor
  
  ;-- Determine OBSERVER:
  observer = 'Name'
  
  ;-- Determine ANGULAR_WIDTH_TIME:
  if exist(timeAngWidthSet) then begin
    timeAngWidthSet=timeAngWidthSet*stepSize/3.
    angular_width_time_day=strarr(2)
    angular_width_time=strarr(2)
    angular_width_file=strarr(2)
    idxDateObs=lonarr(2)
    for i=0,1 do begin
      idxDateObs[i]=(where(timeAngWidthSet[i] le time_sec))[0]
      angular_width_time_day[i]=header_day[idxDateObs[i]]
      angular_width_time[i] = strmid(header_time[idxDateObs[i]],0,8)
      angular_width_file[i] = fileNameList[idxDateObs[i]]
    endfor

  endif else begin
    angular_width_time=['',''] & angular_width_time_day=['',''] & angular_width_file=['','']
  endelse
  
  ;-- ANGULAR_WIDTH_POSI:
  if exist(angWidthSet) then begin
    ang_width_pos=strarr(2)
    ang_width_pos[0]=STRTRIM(string(angWidthSet[0],format='(f6.2)'),1)
    ang_width_pos[1]=STRTRIM(string(angWidthSet[1],format='(f6.2)'),1)
  endif else begin
    ang_width_pos=['','']
  endelse
  
  ;-- Determine ANGULAR_WIDTH:
  if exist(angWidthSet) then begin
    angular_width = STRTRIM(string(abs(angWidthSet[0]-angWidthSet[1]),format='(f6.2)'),1)
  endif else begin
    angular_width = ''
  endelse
  
  ;-- Determine PRE-EVENT_TIME:
  if exist(previousTime) then begin
    previousTime=previousTime*stepSize/3.
    previousTime_day=strarr(1)
    idxDateObs=(where(previousTime le time_sec))[0]
    previousTime_day=header_day[idxDateObs]
    
    previousTime_time=strarr(1)
    previousTime_time = strmid(header_time[idxDateObs],0,8)
  endif else begin
    previousTime_day='' & previousTime_time=''
  endelse
  
  ; PRE-EVENT_FILENAME:
  
  preEventFilename=fileNameList[idxDateObs]
  
  ;-- Set of clicked points
  ;- HEIGHT(SR):
  heightSR=yAxisClickSet/onersun()
  
  ;- DATE:
  date_set=strarr(n_elements(xAxisClickSet[*,0]),n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    for j=0,setnpoints[i]-1 do begin
      idxDateObs=(where(xAxisClickSet[j,i] le time_sec))[0]
      date_set[j,i]=header_day[idxDateObs]
    endfor
  endfor
  
  ;- TIME:
  time_set=strarr(n_elements(xAxisClickSet[*,0]),n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    for j=0,setnpoints[i]-1 do begin
      time_obs_sec = xAxisClickSet[j,i]+time_sec0[0]-long(STRMID(header_day[0],8,2))*60*60*24
      time_obs_hr = floor(time_obs_sec/(60.*60.))
      time_obs_min = floor(time_obs_sec/60.-time_obs_hr*60.)
      sTime_obs_sec = STRTRIM(string(time_obs_sec - time_obs_min*60. - time_obs_hr*60.*60,format='(i2.2)'),1)
      sTime_obs_min = STRTRIM(string(time_obs_min,format='(i2.2)'),1)
      if long(time_obs_hr) ge 24 then time_obs_hr = 0
      sTime_obs_hr = STRTRIM(string(time_obs_hr,format='(i2.2)'),1)
      time_set[j,i] = sTime_obs_hr + ':' + sTime_obs_min + ':' + sTime_obs_sec
    endfor
  endfor
  
  ;- FILENAMELIST:
  filename_set=strarr(n_elements(xAxisClickSet[*,0]),n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    for j=0,setnpoints[i]-1 do begin
      idxDateObs=(where(xAxisClickSet[j,i] le time_sec))[0]
      filename_set[j,i]=fileNameList[idxDateObs]
    endfor
  endfor
  
  ;- 1ST_ORDER_SPEED:
  speedEquation1st=strarr(n_elements(xstart))
  A1=strarr(n_elements(xstart))
  Ea1=strarr(n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    A1[i]=STRTRIM(string(setSpeed1stOrder[i],format='(f8.3)'),1)
    Ea1[i]=STRTRIM(string(setSigma1stOrder[i],format='(f8.3)'),1)
    speedEquation1st[i]= A1[i] +' +- '+ Ea1[i]
  endfor
  
  ;COL(PX)
  
  col_set=lonarr(n_elements(xAxisClickSet[*,0]),n_elements(xstart))
  for i=0,n_elements(xstart)-1 do begin
    for j=0,setnpoints[i]-1 do begin
      col_set[j,i]=(where(xclickpx[j,i] le xaxis_position))[0]
    endfor
  endfor
  
  
  cd,path
  openw,fid,filename
  printf,fid,'#MVI_FILE: ' + mviFileName
  printf,fid,'#PLOT_START_TIME: ' + plot_start_time
  printf,fid,'#FILENAME_START_TIME: ' + fileNameList[0]
  printf,fid,'#PLOT_END_TIME: ' + plot_end_time
  printf,fid,'#FILENAME_END_TIME: ' + fileNameList[n_elements(fileNameList)-1]
  printf,fid,'#DETECTOR: ' + detector
  printf,fid,'#CDELT1: ' + STRTRIM(string(cdelt1,format='(f6.2)'),1)
  printf,fid,'#RSUN: ' + STRTRIM(string(rsun,format='(f6.2)'),1)
  printf,fid,'#N_FILES: ' + STRTRIM(string(n_images),1)
  printf,fid,'#OBSERVER: ' + observer
  printf,fid,'#IMAGE_TYPE: ' + image_type
  printf,fid,'#PRE-EVENT_TIME: ' + previousTime_day + ' ' + previousTime_time
  printf,fid,'#PRE-EVENT_FILENAME: ' + preEventFilename
  printf,fid,'#POSITION_ANGLE: ' + STRTRIM(string(curAngle,format='(f6.2)'),1)
  printf,fid,'#ANGULAR_WIDTH_POSI1: ' + ang_width_pos[0]
  printf,fid,'#ANGULAR_WIDTH_TIME1: ' + angular_width_time_day[0] + ' ' + angular_width_time[0]
  printf,fid,'#ANGULAR_WIDTH_FILE1: ' + angular_width_file[0]
  printf,fid,'#ANGULAR_WIDTH_POSI2: ' + ang_width_pos[1]
  printf,fid,'#ANGULAR_WIDTH_TIME2: ' + angular_width_time_day[1] + ' ' + angular_width_time[1]
  printf,fid,'#ANGULAR_WIDTH_FILE2: ' + angular_width_file[1]
  printf,fid,'#ANGULAR_WIDTH: ' + angular_width
  printf,fid,'#DIST_SUN: ' + STRTRIM(string(long(inRadius)),1)
  printf,fid,'#ATTACHED_IMAGES_1: ' + fnameout1
  printf,fid,'#ATTACHED_IMAGES_2: ' + fnameout4
  printf,fid,'#ATTACHED_IMAGES_3: ' + fnameout5
  printf,fid,'#ATTACHED_IMAGES_4: ' + fnameout6
  printf,fid,'#ATTACHED_IMAGES_5: ' + fnameout7
  printf,fid,'#COMMENT: '
  for i=0,n_elements(xstart)-1 do begin
    printf,fid,'#'
    printf,fid,'#CURVE_N: ' + STRTRIM(string(i+1),1)
    printf,fid,'#DATE-OBS: ' + date_obs[i]
    printf,fid,'#TIME-OBS: ' + time_obs[i]
    printf,fid,'#1ST_ORDER_SPEED: ' + 'V1st(km/s) = ' + speedEquation1st[i]
    printf,fid,'#1ST_ORDER_SPEED: ' + 'V1st(km/s) = A1 +- Ea1'
    printf,fid,'#A1: ' + STRTRIM(string(A1[i],format='(f8.3)'),1)
    printf,fid,'#Ea1: ' + STRTRIM(string(Ea1[i],format='(f8.3)'),1)
    printf,fid,'#2ND_ORDER_SPEED: ' + 'V2nd(km/s) = (' + speedEquation2nd[i]
    printf,fid,'#2ND_ORDER_SPEED: ' + 'V2nd(km/s) = (A2 +- Ea2)*t(s) + (B2 +- Eb2)'
    printf,fid,'#A2: ' + STRTRIM(string(2*A[i],format='(f8.5)'),1)
    printf,fid,'#Ea2: ' + STRTRIM(string(2*Ea[i],format='(f8.5)'),1)
    printf,fid,'#B2: ' + STRTRIM(string(B[i],format='(f9.3)'),1)
    printf,fid,'#Eb2: ' + STRTRIM(string(Eb[i],format='(f9.3)'),1)
    printf,fid,'#N_POINTS: ' + STRTRIM(string(setnpoints[i]),1)
    printf,fid,'#REMARK: ' + remarkSet[i]
    printf,fid,'# HEIGHT(1M Km) HEIGHT(SR)   DATE        TIME    COL(px)  ROW(px)           FILENAME'
    for j = 0,setnpoints[i]-1 do begin
      printf,fid,format='(3X, f6.2, 7X, f6.2, 4X, A10, 3X, A8, 3X, i3, 7X, i3, 4X, A)', $
        yAxisClickSet[j,i]/10.^6.,heightSR[j,i],header_day[col_set[j,i]],strmid(header_time[col_set[j,i]],0,8),round(xclickpx[j,i]),round(yclickpx[j,i]),fileNameList[col_set[j,i]]
    endfor
  endfor
  close,fid
end