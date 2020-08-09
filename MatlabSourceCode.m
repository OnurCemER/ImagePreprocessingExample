function varargout = odev(varargin)
% ODEV MATLAB code for odev.fig
%      ODEV, by itself, creates a new ODEV or raises the existing
%      singleton*.
%
%      H = ODEV returns the handle to a new ODEV or the handle to
%      the existing singleton*.
%
%      ODEV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ODEV.M with the given input arguments.
%
%      ODEV('Property','Value',...) creates a new ODEV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before odev_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to odev_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help odev

% Last Modified by GUIDE v2.5 23-Jun-2020 14:38:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @odev_OpeningFcn, ...
                   'gui_OutputFcn',  @odev_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before odev is made visible.
function odev_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to odev (see VARARGIN)

% Choose default command line output for odev
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes odev wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = odev_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnResimKaydet.
function btnResimKaydet_Callback(hObject, eventdata, handles)
% hObject    handle to btnResimKaydet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[res,path]=uiputfile({'*.jpg';'*.gif';'*.png';'*.ppm';'*.pgm'},'Resmi Kaydet');
f=getframe(handles.axes2);
x=frame2im(f);
imwrite(x,fullfile(path,res));
myicon = x;
msgbox('Resim Kaydedildi','Ýþlem Baþarýlý','custom',myicon);

% --- Executes on button press in btnResimYukle.
function btnResimYukle_Callback(hObject, eventdata, handles)
% hObject    handle to btnResimYukle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[resim,path]=uigetfile({'*.jpg';'*.gif';'*.png';'*.ppm';'*.pgm'}, 'Görüntü Dosyasýný Seç:','Multiselect','On');
oku = fullfile(path,resim);
resim=imread(oku);

axes(handles.axes1);
imshow(resim);
setappdata(0,'resim',resim);
[a,b,c] = size(resim);
sonuc = a*b*c;
kb = sonuc/1024
set(handles.lblOncekiBoyut, 'String', kb);

% --- Executes on button press in btnReset.
function btnReset_Callback(hObject, eventdata, handles)
% hObject    handle to btnReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnMean.
function btnMean_Callback(hObject, eventdata, handles)
% hObject    handle to btnMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

res = getappdata(0,'resim');
b = padarray(res,[1,1]);
b = double(b);
cekirdek = [1,1,1;1,1,1;1,1,1]/9;
output = zeros(size(res));
for i=1:size(b,1)-2
    for j=1:size(b,2)-2
        temp=b(i:i+2,j:j+2).*cekirdek;
        output(i,j)=sum(temp(:));
    end
end
axes(handles.axes2);
imshow(uint8(output));


% --- Executes on button press in btnSobel.
function btnSobel_Callback(hObject, eventdata, handles)
% hObject    handle to btnSobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
img = rgb2gray(res);
db = double(img);
for i=1:size(db,1)-2
    for j=1:size(db,2)-2
        gx=((2*db(i+2,j+1)+db(i+2,j)+db(i+2,j+2))-(2*db(i,j+1)+db(i,j)+db(i,j+2)));
        gy=((2*db(i+1,j+2)+db(i,j+2)+db(i+2,j+2))-(2*db(i+1,j)+db(i,j)+db(i+2,j)));
        img(i,j) = sqrt(gx.^2+gy.^2);
    end
end
axes(handles.axes2);
imshow(uint8(img));
        
% --- Executes on button press in btnGauss.
function btnGauss_Callback(hObject, eventdata, handles)
res = getappdata(0,'resim');
I = double(res);
sigma = 1.75;

sz = 5;
[x,y]=meshgrid(-sz:sz,-sz:sz);

M = size(x,1)-1;
N = size(y,1)-1;
Exp_comp = -(x.^2+y.^2)/(2*sigma*sigma);
Kernel= exp(Exp_comp)/(2*pi*sigma*sigma);

Output=zeros(size(I));

I = padarray(I,[sz sz]);


for i = 1:size(I,1)-M
    for j =1:size(I,2)-N
        Temp = I(i:i+M,j:j+M).*Kernel;
        Output(i,j)=sum(Temp(:));
    end
end
axes(handles.axes2);
imshow(uint8(Output));


% --- Executes on button press in btnMedian.
function btnMedian_Callback(hObject, eventdata, handles)
% hObject    handle to btnMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
matris = ones(3);
output = zeros(size(res)-2);
for i=1:size(res,1)-2
    for j=1:size(res,2)-2
        temp = double(res(i:i+2,j:j+2)).*matris;
        output(i,j)= median(temp(:));
    end
end
axes(handles.axes2);
imshow(uint8(output));

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
res = double(res);
matris = [-1,-1,-1;-1,8,-1;-1,-1,-1];
a = padarray(res,[1 1]);
output=zeros(size(res));
for i=1:size(a,1)-2
    for j=1:size(a,2)-2
        temp = a(i:i+2,j:j+2).*matris;
        output(i,j) = sum(temp(:));
    end
end
axes(handles.axes2);
imshow(uint8(output));

% --- Executes on button press in btnPrewitt.
function btnPrewitt_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res= getappdata(0,'resim'); 
res= rgb2gray(res); 
res= double(res); 
img = zeros(size(res));   
Mx = [-1,0,1; -1,0,1; -1,0,1]; 
My = [-1,-1,-1; 0,0,0; 1,1,1]; 
  
for i = 1:size(res, 1) - 2 
    for j = 1:size(res, 2) - 2 
        Gx = sum(sum(Mx.*res(i:i+2, j:j+2))); 
        Gy = sum(sum(My.*res(i:i+2, j:j+2))); 
        img(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
    end
end
esik = 100;
out = max(img, esik); 
out(out == round(esik)) = 0; 
axes(handles.axes2);
imshow(uint8(out));


% --- Executes on button press in btnRoberts.
function btnRoberts_Callback(hObject, eventdata, handles)
% hObject    handle to btnRoberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
res = rgb2gray(res);
out =edge(res,'Roberts',0.08);
axes(handles.axes2)
imshow(out);

% --- Executes on button press in btnAsinma.
function btnAsinma_Callback(hObject, eventdata, handles)
% hObject    handle to btnAsinma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
res = im2bw(res);
matris = [0,1,0;1,1,1;0,1,0];
satir = padarray(res,[0 1],1);
q = false(size(res));
for i=1:size(satir,1)-2
    for j=1:size(satir,2)-2
        L = satir(i:i+2,j:j+2);
        K = find(matris==1);
        if(L(K)==1)
            q(i,j)=1;
        end
    end
end
axes(handles.axes2);
imshow(q);

% --- Executes on button press in btnGenisleme.
function btnGenisleme_Callback(hObject, eventdata, handles)
% hObject    handle to btnGenisleme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
res = im2bw(res);
matris = [0,1,0;1,1,1;0,1,0];
c= padarray(res,[0 3]);
q = false(size(res));
for i=1:size(c,1)-2
    for j=1:size(c,2)-2
        q(i,j) = sum(sum(matris&c(i:i+2,j:j+2)));
    end
end
axes(handles.axes2);
imshow(q);

% --- Executes on button press in btnKapama.
function btnKapama_Callback(hObject, eventdata, handles)
res = getappdata(0,'resim');
res = im2bw(res);

c= padarray(res,[0 3]);
matris = [1,1,1;1,1,1;1,1,1];
q = false(size(res));

for i=1:size(c,1)-2
    for j=1:size(c,2)-2
        q(i,j) = sum(sum(matris&c(i:i+2,j:j+2)));
    end
end
matris = [1,1,0];
satir = padarray(res,[0 1],1);
q = false(size(res));
for i=1:size(satir,1)-2
    for j=1:size(satir,2)-2
        z = satir(i:i+2,j:j+2);
        w = find(matris==1);
        if(z(w)==1)
            q(i,j)=1;
        end
    end
end
axes(handles.axes2);
imshow(q);

% --- Executes on button press in btnAcma.
function btnAcma_Callback(hObject, eventdata, handles)
res = getappdata(0,'resim');
res = im2bw(res);
matris = [1,1,0];
satir = padarray(res,[0 1],1);
q = false(size(res));

for i=1:size(satir,1)-2
    for j=1:size(satir,2)-2
        z = satir(i:i+2,j:j+2);
        w = find(matris==1);
        if(z(w)==1)
            q(i,j)=1;
        end
    end
end
c= padarray(res,[0 3]);
matris = [1,1,1;1,1,1;1,1,1];
q = false(size(res));
for i=1:size(c,1)-2
    for j=1:size(c,2)-2
        q(i,j) = sum(sum(matris&c(i:i+2,j:j+2)));
    end
end
axes(handles.axes2);
imshow(q);
% --- Executes on button press in btnEsikleme.
function btnEsikleme_Callback(hObject, eventdata, handles)
% hObject    handle to btnEsikleme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
val = str2double(get(handles.txtEsik,'String'));
if isempty(handles.txtEsik.String)
    f = msgbox("Eþik Deðeri Boþ Býrakýlamaz.");
    uiwait(warndlg(f));
elseif isnan(val)
   f = msgbox("Eþik Deðeri Bir Sayý Olmalýdýr. ÖRN: 0.5 || 50");
   uiwait(warndlg(f));
elseif val<0
   f = msgbox("Eþik Deðeri Negatif Bir Sayý Olamaz.");
   uiwait(warndlg(f));
end
try
    imge=rgb2gray(res);
catch
    imge = im2double(res);
end
axes(handles.axes2);
imshow(imge>val);


% --- Executes on button press in btnLog.
function btnLog_Callback(hObject, eventdata, handles)
% hObject    handle to btnLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
doubleRes = im2double(res);
img = doubleRes;
[r,c] = size(doubleRes);
val = str2double(get(handles.txtLog,'String'));
if isempty(handles.txtLog.String) 
   f = msgbox("Katsayý Deðeri Boþ Býrakýlamaz.");
   uiwait(warndlg(f));
elseif isnan(val)
    f = msgbox("Katsayý Bir Sayý Olmalýdýr. ÖRN: 5");
    uiwait(warndlg(f));
elseif(val<=0)
    f = msgbox("Katsayý 0 veya daha küçük olamaz. ÖRN: 5");
    uiwait(warndlg(f));
end
for i=1:r
    for j=1:c
        img(i,j)=val*log(1+doubleRes(i,j));
    end
end
axes(handles.axes2);
imshow(img);

% --- Executes on button press in btnKontrastGerme.
function btnKontrastGerme_Callback(hObject, eventdata, handles)
% hObject    handle to btnKontrastGerme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
maxRes=max(max(res));
minVal = str2double(get(handles.txtMinimum,'String'));
maxVal = str2double(get(handles.txtMaximum,'String'));
if isempty(handles.txtMinimum.String) || isempty(handles.txtMaximum.String) 
   f = msgbox("Deðerlerin ikisi de doldurulmalýdýr.");
   uiwait(warndlg(f));
elseif isnan(minVal) || isnan(maxVal)
    f = msgbox("Giriþler Bir Sayý Olmalýdýr.");
    uiwait(warndlg(f));
elseif(minVal<=0)
    f = msgbox("Minimum sayý 0 veya daha küçük olamaz.");
    uiwait(warndlg(f));
end
maxRes=maxVal ;
maxRes=double(maxRes);
minRes = min(min(res));
minRes = minVal;
minRes = double(minRes);
o = maxRes/255; o2 = minRes/255;
out = imadjust(res,[o2 o],[0 1]);
axes(handles.axes2);
imshow(out);

% --- Executes on button press in btnGriSeviye.
function btnGriSeviye_Callback(hObject, eventdata, handles)
% hObject    handle to btnGriSeviye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim');
try
  [satir, sutun, renkSayisi] = size(res);
  if renkSayisi  == 3
      red = res(:, :, 1);
      green = res(:, :, 2);
      blue = res(:, :, 3);
      griSeviye = .299*double(red) + ...
                  .587*double(green) + ...
                  .114*double(blue);
      griSeviye = uint8(griSeviye);
      axes(handles.axes2);
      imshow(griSeviye);
  else
      griSeviye = res;
      axes(handles.axes2);
      imshow(griSeviye);
  end
catch ME
  f = msgbox("Bir Hata Oluþtu!");
    uiwait(warndlg(f));
end

function txtEsik_Callback(hObject, eventdata, handles)
% hObject    handle to txtEsik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEsik as text
%        str2double(get(hObject,'String')) returns contents of txtEsik as a double


% --- Executes during object creation, after setting all properties.
function txtEsik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEsik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res =  getappdata(0,'resim');
res = rgb2gray(res)
max_satir = size(res,1);
max_sutun = size(res,2);
histogram = zeros([1 256]);
cumulative_histogram = zeros([1 256]);
for satir=1:max_satir
   for sutun=1:max_sutun
    for sayac=1:256
      if(res(satir,sutun) == sayac-1)
         histogram(sayac) = histogram(sayac) + 1;
          break;
         end
      end
   end 
end
current_value = 0;
for sayac=1:256
    current_value = current_value + histogram(sayac);
    cumulative_histogram(sayac) = current_value;
end
normal_hist = zeros([1 256]);
cdf_min = min(cumulative_histogram);
for sayac=1:256 
normal_hist(sayac) = cumulative_histogram(sayac) - cdf_min;
normal_hist(sayac) = normal_hist(sayac) / ((max_satir*max_sutun)-cdf_min);
normal_hist(sayac) = round(normal_hist(sayac) * 255);
end
esitleme_img = zeros([max_satir max_sutun]);
for satir=1:max_satir
   for sutun=1:max_sutun
      for sayac=1:256
      if(res(satir,sutun) == (sayac-1))
        esitleme_img(satir,sutun) = normal_hist(sayac);
      break;
       end
     end
   end
end

axes(handles.axes2);
imshow(uint8(esitleme_img));


function txtLog_Callback(hObject, eventdata, handles)
% hObject    handle to txtLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLog as text
%str2double(get(hObject,'String')) returns contents of txtLog as a double


% --- Executes during object creation, after setting all properties.
function txtLog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnHuffman.
function btnHuffman_Callback(hObject, eventdata, handles)
% hObject    handle to btnHuffman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res = getappdata(0,'resim')
I = rgb2gray(res);
[m,n]=size(I);
toplam=m*n;
cnt=1;
sigma=0;

for i=0:255
    k=I==i;
    say(cnt)=sum(k(:))
    tk(cnt)=say(cnt)/toplam;
    sigma=sigma+tk(cnt);
    cumtk(cnt)=sigma;
    cnt=cnt+1;
end

symbols = [0:255];
dict = huffmandict(symbols,tk);
vec_size = 1;
for p = 1:m
    for q = 1:n
        yeniVek(vec_size) = I(p,q);
        vec_size = vec_size+1;
    end
end
hufS = huffmanenco(yeniVek,dict);
hufD = huffmandeco(hufS,dict);
decHufman = uint8(hufD);
satir=sqrt(length(decHufman));
decOku=satir;
satir1 = 1;
sutun1 = 1;
vektor = 1;
for x = 1:m
    for y = 1:n
    sonHal(x,y)=decHufman(vektor);
    sutun1 = sutun1+1;
    vektor = vektor + 1;
    end
    satir1 = satir1+1;
end
[deco, map] = gray2ind(sonHal,256);
out = ind2rgb(deco,map);
axes(handles.axes2);
imshow(out);
[gri,res,c] = size(out);
sonuc = gri*res*c;
kb = sonuc/1024;
set(handles.lblSonrakiBoyut, 'String', sonuc);

% --- Executes on button press in btnDCT.
function btnDCT_Callback(hObject, eventdata, handles)
% hObject    handle to btnDCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resim = getappdata(0,'resim');
img =resim;
img=double(img);
satirSay=size(img,1);
kolonSay=size(img,2);
strBlok=floor(kolonSay/8);
klnBlok=floor(satirSay/8);
mtr1=[];
for k=1:3 
    for i=1:satirSay        
     satir_DCT=dct(img(i,:,k));
     mtr1(i,:,k)=idct(satir_DCT(1:strBlok),satirSay); 
   end
end

mtr2=[];
for k=1:3    
   for i=1:kolonSay       
    kolon_DCT=dct(mtr1(:,i,k));
    mtr2(:,i,k)=idct(kolon_DCT(1:klnBlok),kolonSay); 
  end
end

axes(handles.axes2);
imshow(uint8(mtr2));
[a,b,c] = size(mtr2);
sonuc = a*b*c;
kb = sonuc/1224;
set(handles.lblSonrakiBoyut, 'String', kb);

% --- Executes on button press in btnAritmetik.
function btnAritmetik_Callback(hObject, eventdata, handles)
% hObject    handle to btnAritmetik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnVektor.
function btnVektor_Callback(hObject, eventdata, handles)
% hObject    handle to btnVektor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
res=getappdata(0,'resim');
cbs=256;
bs=4;
[h,w]=size(res);
res=double(res);
cb=zeros(bs^2,cbs);
ortak=im2col(res,[bs bs],'sliding');
[m,n]=size(ortak);
q=randperm(n);
cb(:)=ortak(:,q(1:cbs));

index=zeros(1,n);
for t =1:30
    % temp=zeros(size(cb));
    for i = 1:n
        temp=repmat(ortak(:,i),1,cbs);
        [val,ind]=min(sum((temp-cb).^2));
        index(i)=ind;
    end
    q=1:n;
    for i =1:cbs
        cb(:,i)=round(mean(ortak(:,q(index==i))')');    
    end
    t
end
save cb cb
load cb
res=double(res);
ortak=im2col(res,[bs bs],'distinct');
[m,n]=size(res);
index=zeros(m/bs,n/bs);
for i = 1:numel(index)
    temp=repmat(ortak(:,i),1,cbs);
    [val,ind]=min(sum((temp-cb).^2));
    index(i)=ind;        
end
ortak=cb(:,index(:));
out=uint8(col2im(ortak,[bs bs],[h,w],'distinct'));
axes(handles.axes2);
imshow(out);
[a,b] = size(out);
sonuc = a*b;
kb = sonuc/2048;
set(handles.lblSonrakiBoyut, 'String', kb);
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMinimum_Callback(hObject, eventdata, handles)
% hObject    handle to txtMinimum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMinimum as text
%        str2double(get(hObject,'String')) returns contents of txtMinimum as a double


% --- Executes during object creation, after setting all properties.
function txtMinimum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMinimum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMaximum_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaximum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaximum as text
%        str2double(get(hObject,'String')) returns contents of txtMaximum as a double


% --- Executes during object creation, after setting all properties.
function txtMaximum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaximum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
