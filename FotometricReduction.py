'''
----------------------------------TR----------------------------------
Bu kod, NY Vir star'ın fotometrik indirgemesi yapılıp ışık eğrisi hesaplanması için kullanılmıştır. İşlemler şöyle yapılmıştır:

1. NY Vir görüntüleri için veri indirgemesi yapıldı (master bias ve master flat).
2.NY Vir görüntüleri için ekstra görüntü temizleme ve trim yapıldı.
3. Trimlenmiş görüntüler hizalandı.
4. NY Vir ve ona hem parlaklık hem de konum olarak yakın bir yıldızın piksel konumu seçilerek fotometri yapıldı.
5. Fotometri yapılan NY Vir'in fark ışık ölçümü yapılarak ışık eğrisi çizildi ve minimum zamanı belirlenmesi için uygun bir fit modeli çizildi.
6. Işık eğrisi görüntüsü ve çekilen fit modelinin sonuçları kaydedildi.

Bu kod, veri indirgeme, görüntü işleme ve fotometri analizi yaparak, astronomi araştırmalarında kullanılan standart bir yöntem olan "light curve" analizini gerçekleştirir.
Işık eğrisi, gözlemlenen bir yıldızın ışık şiddetinin zamana göre değişimini gösteren bir grafiktir.
Bu analiz, gözlemlenen yıldızın parlaklık değişimlerinin incelenmesi ve yıldızın fiziksel özellikleri hakkında bilgi edinilmesi için önemlidir.


----------------------------------ENG----------------------------------
This code was used for the photometric reduction of an observed star and calculation of its light curve.
The sequence of operations performed is as follows (The observed star used as an example is NY Vir star.):

1. Data reduction was performed for NY Vir images (Master Bias and Master Flat).
2. Additional image cleaning and trimming were performed for NY Vir images.
3. Trimmed images were aligned.
4. Photometry was performed by selecting the pixel position of NY Vir and a nearby star that is similar in brightness and position.
5. The light curve of NY Vir was drawn by measuring its differential magnitude in the Fark Işık band, and an appropriate fit model was drawn to determine the minimum time.
6. The results of the light curve image and the drawn fit model were recorded.
'''
# Gerekli kütüphaneler / Required libraries
import os
from astropy.nddata.ccddata import CCDData
from ccdproc.combiner import combine
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from astropy.io import fits
import ccdproc
from astropy import units as u
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter, MultipleLocator
from lmfit.models import GaussianModel, LorentzianModel, VoigtModel, PolynomialModel, ConstantModel, ExponentialGaussianModel, ExponentialModel
import scipy.stats as stats
from sklearn.metrics import r2_score
import pandas as pd
import astroalign as aa
from astropy.io.fits.convenience import update
from astropy.nddata import CCDData
from numpy.lib.twodim_base import tri
from ccdproc.core import ccd_process
from astropy.visualization import (
    PowerStretch, ZScaleInterval, ImageNormalize)
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import aperture_photometry
from astropy.visualization import simple_norm
from astropy.stats import sigma_clipped_stats
from photutils.utils import calc_total_error
from photutils.centroids import centroid_sources, centroid_com, centroid_2dg
import matplotlib.patches as pat


# FITS görüntüleri klasörleri (Bu bölümde çekilecek görüntülerin tipine göre klasör oluşturulmuştur.):
# FITS image directories (Directories were created according to the type of images to be taken in this section).
sci = []  # SCI (Light) Dosyaları / SCI (Light) Files
bias = []  # Bias Dosyaları / Bias Files
flat = []  # Flat Dosyaları / Flat Files
bias_ss = []  # Bias Standart Sapma / Bias Standard Deviation
flat_ort = []  # Flat Ortalama / Flat Mean
m_sub_flat = []  # Bias Çıkartılmış Flatlerin Listesi / List of Bias Subtracted Flats
calibrated_sci = []  # Kalibre edilmiş SCI görüntüleri / Calibrated SCI Images
x = 0


def fileManagement(sci, bias, flat):

    fts = glob('C:\\path\\NYVir\\*.fits')
    fts.extend(glob('C:\\path\\NYVIR\\BDF\\*.fits'))
    for f in fts:
        hdu = fits.open(f, mode='update')
        hdr = hdu[0].header
        ityp = hdr['IMAGETYP']

        if ityp == 'Light':
            sci.append(f)
        elif ityp == 'Bias':
            bias_std = hdu[0].data
            std_bias = np.std(bias_std)
            bias_ss.append(std_bias)
            if 10 < std_bias < 20:
                bias.append(CCDData(data=hdu[0].data, meta=hdr, unit="adu"))
        elif ityp == 'Flat':
            flat_mean = hdu[0].data
            mean_flat = np.mean(flat_mean)
            flat_ort.append(mean_flat)
            if 10000 < mean_flat < 40000:
                flat.append(CCDData(data=hdu[0].data, meta=hdr, unit="adu"))
        else:
            print('There is a different image type: {}, {}'.format(f, ityp))
        if hdr['XBINNING'] != 2:
            print('There is a different binning: {}, {}'.format(
                f, hdr['XBINNING']))


fileManagement(sci, bias, flat)

'''
----------------------------------TR----------------------------------
Bu kod, belirtilen yoldaki tüm .fits dosyalarını seçer ve dosya başlıklarındaki IMAGETYP etiketine bakarak dosyaları uygun klasörlere yerleştirir.
sci, bias ve flat adında üç farklı liste oluşturulur. 'sci' listesi, IMAGETYP = 'Light' olan .fits dosyalarını içerir.
'bias' listesi, IMAGETYP = 'Bias' olan .fits dosyalarını ve her bias dosyasının standart sapmasını içerir.
'flat' listesi, IMAGETYP = 'Flat' olan .fits dosyalarını ve her flat dosyasının ortalamasını içerir.
Dosyaların farklı bir IMAGETYP veya binning değeri varsa, bu dosyaların konumu ekrana yazdırılır.
Bu işlemlerden sonra, 'fileManagement' fonksiyonu ile SCI, BIAS ve FLAT dosyaları ayrı ayrı düzenlenir ve bilimsel fotometrik analizler için hazır hale getirilir.
Komutların açıklamaları:
'glob' : Verilen bir kalıp eşleşmesi için tüm dosyaları seçmek için kullanılan bir fonksiyondur.
'fits.open' : .fits dosyalarını açmak için kullanılan bir fonksiyondur.
'hdu' : Açılan .fits dosyasındaki tüm verileri ve başlıkları içeren bir liste nesnesidir.
'hdr' : Her bir dosyanın header bilgilerini açmak için kullanılan bir değişkendir.
'ityp' : Header bilgilerindeki IMAGETYP bölümünü tutan bir değişkendir.
'CCDData' : Bir CCD kamera tarafından alınan verileri saklamak için kullanılan bir sınıftır.
'np.std' : Standart sapmayı hesaplamak için kullanılan bir fonksiyondur.
'np.mean' : Ortalama değeri hesaplamak için kullanılan bir fonksiyondur.
'print' : Belirtilen mesajı ekrana yazdırmak için kullanılan bir fonksiyondur.
'if-elif-else' : Koşul ifadeleri için kullanılan yapıdır.
'fileManagement' : Bilimsel fotometrik analizler için SCI, BIAS ve FLAT dosyalarını hazırlamak için kullanılan bir fonksiyondur.

----------------------------------ENG----------------------------------
This code selects all .fits files in the specified directory and places them in appropriate folders based on the IMAGETYP tag in the file headers.
Three different lists are created: sci, bias, and flat. The 'sci' list contains .fits files with IMAGETYP='Light'.
The 'bias' list contains .fits files with IMAGETYP='Bias' and the standard deviation of each bias file.
The 'flat' list contains .fits files with IMAGETYP='Flat' and the average of each flat file.
If the files have a different IMAGETYP or binning value, their location is printed to the screen.
After these operations, the 'fileManagement' function organizes the SCI, BIAS, and FLAT files separately and prepares them for scientific photometric analysis.
Descriptions of the commands:
'glob': A function used to select all files for a given pattern match.
'fits.open': A function used to open .fits files.
'hdu': A list object containing all data and headers in the opened .fits file.
'hdr': A variable used to open the header information for each file.
'ityp': A variable holding the IMAGETYP section in the header information.
'CCDData': A class used to store data taken by a CCD camera.
'np.std': A function used to calculate standard deviation.
'np.mean': A function used to calculate mean.
'print': A function used to print a specified message to the screen.
'if-elif-else': A structure used for conditional statements.
'fileManagement': A function used to prepare SCI, BIAS, and FLAT files for scientific photometric analysis.
'''


def master_combine(bias, flat, m_sub_flat):
    if os.path.isfile('master_bias.fits'):
        os.remove('master_bias.fits')
    mstr_bias = combine(bias, method='median', unit="adu")
    mstr_bias.write('D:\\path\\master_bias.fits',
                    overwrite=True, output_verify='silentfix')

    for flt in flat:
        sub_flat = ccdproc.subtract_bias(flt, mstr_bias)
        m_sub_flat.append(sub_flat)
    if os.path.isfile('master_flat.fits'):
        os.remove('master_flat.fits')
    mstr_flat = combine(m_sub_flat, method='average', unit="adu")
    mstr_flat.write('C:\\path\\master_flat.fits',
                    overwrite=True, output_verify='silentfix')


master_combine(bias, flat, m_sub_flat)


'''
----------------------------------TR----------------------------------
Bu kod, 'bias' ve 'flat' adlı iki liste ve 'm_sub_flat' adlı bir boş liste kullanarak, kalibre edilmiş SCI görüntüleri oluşturur.

'master_combine' adlı bir fonksiyon, 'bias' listesinden master bias görüntüsünü oluşturur.
Daha sonra 'flat' listesindeki her flat görüntüsünden master_bias görüntüsü çıkarılır ve sonuçta 'm_sub_flat' listesinde kalibre edilmiş flat görüntüleri elde edilir.
Ardından, 'm_sub_flat' listesindeki kalibre edilmiş flat görüntüleri, 'combine' fonksiyonu ile kombine edilerek master flat görüntüsü oluşturulur.

Son olarak, her bir SCI görüntüsünden master_bias ve master_flat görüntüleri çıkarılarak, kalibre edilmiş SCI görüntüleri elde edilir.
Bu kalibre edilmiş SCI görüntüleri, 'calibrated_sci' listesine eklenir.

----------------------------------ENG----------------------------------
This code creates calibrated SCI images using two lists named 'bias' and 'flat' and an empty list named 'm_sub_flat'.

A function named 'master_combine' creates a master bias image from the 'bias' list.
Then, the master_bias image is subtracted from each flat image in the 'flat' list to obtain calibrated flat images in the 'm_sub_flat' list.
Next, the calibrated flat images in the 'm_sub_flat' list are combined using the 'combine' function to create a master flat image.

Finally, calibrated SCI images are obtained by subtracting the master_bias and master_flat images from each SCI image.
These calibrated SCI images are appended to the 'calibrated_sci' list.
'''

master_bias = CCDData.read('C:\\path\\master_bias.fits', unit='adu')  

master_flat = CCDData.read('C:\\path\\master_flat.fits', unit='adu')


def calibrate(sci):
    for _sci in sci:
        if 'bf_' in _sci:
            continue
    ccs = CCDData.read(_sci, unit="adu")
    if np.min(ccs) < 0:
        print(_sci, np.min(ccs), np.max(ccs))
    ccs = ccd_process(ccs, error=True,
                      gain=0.57*u.electron/u.adu,
                      readnoise=4.19*u.electron,
                      master_bias=master_bias,
                      master_flat=master_flat,
                      exposure_key='exposure',
                      exposure_unit=u.second,
                      gain_corrected=False)
    cname = '{}/bf_{}'.format(os.path.dirname(_sci), os.path.basename(_sci))
    if not os.path.isfile(cname):
        ccs.write(cname, output_verify='silentfix')


calibrate(sci)

'''
----------------------------------TR----------------------------------
Bu kod, SCI klasöründeki bilimsel görüntüleri kalibre eder.
Fonksiyon, SCI klasöründeki her görüntüyü okur ve 'bf_' ile başlamayanları kalibre eder.
Görüntülerin minimum piksel değerleri kontrol edilir ve eğer negatif bir değer varsa, bu görüntünün adı ve minimum piksel değeri yazdırılır.
Görüntüler, 'ccd_process' işlevi kullanılarak işlenir.
İşlem sırasında, görüntülerin gain ve readnoise bilgileri belirtilir, master bias ve master flat düzeltmeleri uygulanır, pozlama süresi değişkeni tanımlanır ve çıktı birimleri belirlenir.
Kalibre edilmiş görüntüler 'bf_' ön adı ile kaydedilir.

----------------------------------ENG----------------------------------
This code calibrates the scientific images in the SCI directory.
The function reads each image in the SCI directory and calibrates those that do not start with 'bf_'.
The minimum pixel values of the images are checked, and if there is a negative value, the name of the image and the minimum pixel value are printed.
The images are processed using the 'ccd_process' function.
During processing, the gain and readnoise information of the images are specified, master bias and master flat corrections are applied, the exposure time variable is defined, and the output units are determined.
The calibrated images are saved with the prefix 'bf_'.
'''


img_list = glob('C:\\path\\NYVIR\\bf_NYVir_*_Empty.fits')


def trim(img_list):
    for _img in img_list:
        img_hdu = fits.open(_img, mode='update')
        hdr = img_hdu[0].header
        data = img_hdu[0].data
        ccs = CCDData.read(_img, unit="adu")
        if np.min(ccs) < 0:
            print(_img, np.min(ccs), np.max(ccs))
        ccs = ccd_process(ccs, trim='[20:2028, 20:2028]',
                          exposure_key='exposure',
                          exposure_unit=u.second,
                          gain_corrected=0.57*u.electron/u.adu)
        cname = '{}/t{}'.format(os.path.dirname(_img), os.path.basename(_img))
        ccs.write(cname, overwrite=True)


trim(img_list)

'''
----------------------------------TR----------------------------------
Bu kod, 'img_list' adlı bir dosya listesindeki tüm 'bf_NYVir_*_Empty.fits' dosyalarını kırparak yeniden kaydeder.
Fonksiyon, her dosyayı açar ve fits dosyası biçiminde kaydedilir.
'ccd_process' işlevi kullanılarak görüntüler kırpılır ve işlenir.
İşlem sırasında, görüntülerin trim sınırları belirlenir, pozlama süresi değişkeni tanımlanır ve çıktı birimleri belirlenir.
Kırpılmış ve yeniden işlenmiş görüntüler 't' ön adı ile kaydedilir.


----------------------------------ENG----------------------------------
This code trims and saves all 'bf_NYVir_*_Empty.fits' files in a file list named 'img_list'.
The function opens each file and saves it in the fits file format.
The images are trimmed and processed using the 'ccd_process' function.
During processing, the trimming limits of the images are defined, exposure time variables are specified, and output units are defined.
The trimmed and processed images are saved with the prefix 't'.
'''

img_list = glob('D:\\path\\NYVIR\\tbf_NYVir_*_Empty.fits')

target = 'D:\\path\\NYVIR\\tbf_NYVir_0153_Empty.fits'

target_hdu = fits.open(target)
target_data = target_hdu[0].data
target_data = target_data.astype('float32')

for _img in img_list:
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    img_hdu = fits.open(_img, mode='update')
    print(_img)
    img_data = img_hdu[0].data
    img_data = img_data.astype('float32')
    img_aligned, footprint = aa.register(img_data, target_data,
                                         detection_sigma=5, min_area=5)
    img_hdu.data = img_aligned
    cname = '{}/a{}'.format(os.path.dirname(_img), os.path.basename(_img))
    img_hdu.writeto(cname, overwrite=True)
    norm = ImageNormalize(img_aligned, interval=ZScaleInterval(),
                          stretch=PowerStretch(2))
    ax.imshow(img_aligned[1073:1273, 682:882], norm=norm, cmap='Greys')
    circle = plt.Circle((100, 100), 9., color='r', fill=False)

    ax.add_artist(circle)
    fig.savefig("D:\\path\\myfig\\{}.png".format(os.path.basename(_img)))


'''
----------------------------------TR----------------------------------
Bu kod, 'img_list' adlı dosya listesindeki tüm 'tbf_NYVir_*_Empty.fits' dosyalarını hizalar ve kaydeder.
İlk olarak, 'target' adlı bir referans dosyası seçilir ve bu dosyanın verileri okunur.
Daha sonra, her görüntü dosyası açılır ve verileri 'float32' veri tipine dönüştürülür. 'Astroalign' kütüphanesi kullanılarak hizalanırlar ve işlenirler.
İşlem sırasında, algılama sigma değeri ve minimum alanı gibi parametreler belirlenir. Hizalanmış görüntüler, 'a' ön adı ile kaydedilir.
Son olarak, bir kesit görüntüsü oluşturulur ve bir daire çizilir. Bu görüntüler kaydedilir ve işlenmiş veriler analiz edilmek üzere görselleştirilir.

----------------------------------ENG----------------------------------
This code aligns and saves all 'tbf_NYVir_*_Empty.fits' files in the 'img_list' file list.
First, a reference file named 'target' is selected and its data is read. Then, each image file is opened and its data is converted to the 'float32' data type.
They are aligned and processed using the 'Astroalign' library. During the process, parameters such as detection sigma value and minimum area are determined.
Aligned images are saved with the prefix 'a'. Finally, a section image is created and a circle is drawn.
These images are saved, and the processed data is visualized for further analysis.
'''

coords = [(800, 1184), (905, 1298)]

img_list = glob('C:\\path\\NYVIR\\atbf_NYVir_*_Empty.fits')

tb = open('C:\\path\\tm.csv', 'w')
tb.write('JD,exptime,x1,y1,mag1,err1,x2,y2,mag2,err2\n')

for _i, target in enumerate(img_list):
    data = fits.getdata(target)
    exp_time = fits.getheader(target)['EXPTIME']
    JD = fits.getheader(target)['JD']
    xx, yy = centroid_sources(data, np.transpose(coords)[0],
                              np.transpose(coords)[1],
                              box_size=5,
                              centroid_func=centroid_2dg)
    coords = [(_x, _y) for _x, _y in zip(xx, yy)]
    aperture = CircularAnnulus(coords, r_in=10, r_out=15)

    annulus_masks = annulus_aperture.to_mask(method='center')
    bkg_median = []
    sigma_all = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        mean, median_sigclip, std = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
        sigma_all.append(std)
    bkg_median = np.array(bkg_median)
    sigma_all = np.array(sigma_all)
    error = calc_total_error(data, np.mean(sigma_all), exp_time)

    phot = aperture_photometry(data, aperture, error=error)
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * aperture.area
    phot['aper_sum_bkgsub'] = phot['apertre_sum'] - phot['aper_bkg']
    phot['mag'] = 25.-2.5*np.log10(phot['aper_sum_bkgsub'].data / exp_time)
    phot['mag_err'] = phot['aperture_sum_err'] / \
        phot['aper_sum_bkgsub'] * phot['mag']
    tb.write('{},{},{},{},{},{},{},{},{},{}\n'.format(
        JD, exp_time,
        phot['xcenter'][0], phot['ycenter'][0],
        phot['mag'][0], phot['mag_err'][0],
        phot['xcenter'][1], phot['ycenter'][1],
        phot['mag'][1], phot['mag_err'][1]))
    fig, ax=plt.subplots(1, 1, figsize=(18, 18))
    norm=ImageNormalize(data[1100:1400, 700:1050], interval=ZScaleInterval(),
                          stretch=PowerStretch(2))
    ax.imshow(data[1100:1400, 700:1050], norm=norm, cmap='Greys')
    circle=plt.Circle((xx[0]-700, yy[0]-1100), 5., color='r', fill=False)
    ax.add_artist(circle)
    circle=plt.Circle((xx[1]-700, yy[1]-1100), 5., color='r', fill=False)
    ax.add_artist(circle)
    circle=plt.Circle((xx[0]-700, yy[0]-1100), 5., color='r', fill=False)
    ax.add_artist(circle)
    circle=plt.Circle((xx[1]-700, yy[1]-1100), 5., color='r', fill=False)
    ax.add_artist(circle)
    circle=plt.Circle((xx[0]-700, yy[0]-1100), 10., color='g', fill=False)
    ax.add_artist( circle )
    circle = plt.Circle((xx[1]-700, yy[1]-1100 ),10. ,color='g',fill = False )
    ax.add_artist( circle )
    circle = plt.Circle((xx[0]-700, yy[0]-1100 ),15. ,color='g',fill = False )
    ax.add_artist( circle )
    circle = plt.Circle((xx[1]-700, yy[1]-1100 ),15. ,color='g',fill = False )
    ax.add_artist( circle )
    fname=os.path.basename(target)
    plt.savefig('C:\\path\\png\\{}'.format(fname.replace('fits',''))) 

tb.close()


'''
----------------------------------TR----------------------------------
Bu kod, astronomi alanında kullanılan bir işlem olan astrofotometri işlemini gerçekleştiriyor. 
Astrofotometri, yıldızların, gezegenlerin veya diğer astronomik nesnelerin ışık yoğunluklarını ölçen ve analiz eden bir tekniktir.

Kod, FITS formatındaki görüntü dosyalarını işliyor. 
İlk olarak, glob fonksiyonu kullanılarak belirtilen klasördeki uygun dosyaların listesi oluşturuluyor. 
Daha sonra, bu dosyaların içindeki veriler fits.getdata() ve fits.getheader() fonksiyonları kullanılarak çekiliyor. 
Burada EXPTIME ve JD (Julian Tarih) gibi veriler başlık bilgilerinden alınıyor.

Daha sonra, veriler üzerinde işlem yapmak için bir for döngüsü başlıyor. 
Her iterasyonda, belirli bir görüntü dosyası işleniyor. 
centroid_sources fonksiyonu, veriler içindeki yıldızların merkezlerini bulmak için kullanılıyor. 
Burada, box_size parametresi yıldız merkezlerini belirlemek için kullanılan kutunun boyutunu belirtiyor. 
centroid_2dg fonksiyonu, 2 boyutlu Gauss profili kullanarak yıldızların merkezlerini bulmak için kullanılıyor.

Daha sonra, yıldızların merkezlerini temel alarak bir aperture (yansıtıcılık) oluşturuluyor. 
Bu aperture yıldızların etrafındaki halka şeklindeki alanları kapsayan bir halkadır.
CircularAnnulus fonksiyonu bu halka şeklindeki alanları tanımlamak için kullanılır. 
annulus_aperture.to_mask fonksiyonu, halka alanlarını görüntü verileri üzerinde maskeleme yapmak için kullanılır.

sigma_clipped_stats fonksiyonu, halka alanındaki piksellerin medyan değerini hesaplar. 
Bu değer daha sonra aper_bkg (aper görüntüsü için arka plan) hesaplamak için kullanılır. 
calc_total_error fonksiyonu, hesaplanan arka plan değeri ile birlikte hata değerini hesaplamak için kullanılır.

aperture_photometry fonksiyonu, yıldızların içindeki ışık yoğunluğunu hesaplamak için kullanılır. 
Bu işlem sonucunda aperture_sum, aperture_sum_err gibi değerler hesaplanır. 
Bu değerler kullanılarak yıldızların ışık yoğunluğu (magnitüd) ve hata değerleri hesaplanır. Bu değerler daha sonra .csv dosyasına yazdırılır.

Son olarak, her görüntü için yıldızların merkezlerini ve halka alanlarını gösteren bir matplotlib grafiği oluşturulur.
Daha Sonra "tb.close()" kodu ile dosya kapatılmaktadır. 
Bu işlem, dosyayı yazma işleminden sonra düzgün bir şekilde kapatır ve işlem tamamlandığında açık dosya kaynaklarının serbest bırakılmasını sağlar. 
Bu, programın gereksiz bellek tüketiminden kaçınmasına ve daha düzgün bir şekilde çalışmasına yardımcı olur.


----------------------------------ENG----------------------------------
This code performs an operation used in the field of astronomy called astrophotometry. 
Astrophotometry is a technique that measures and analyzes the light intensities of stars, planets, or other astronomical objects.

The code processes image files in FITS format. 
First, the glob function is used to create a list of suitable files in the specified folder. 
Then, the data within these files is extracted using the fits.getdata() and fits.getheader() functions. 
Here, data such as EXPTIME and JD (Julian Date) are obtained from the header information.

Next, a for loop begins to process the data. In each iteration, a specific image file is processed. 
The centroid_sources function is used to find the centers of the stars within the data. 
Here, the box_size parameter specifies the size of the box used to determine the star centers. 
The centroid_2dg function is used to find the star centers using a 2-dimensional Gaussian profile.

Next, an aperture is created based on the star centers. This aperture is a circular annulus that covers the areas around the stars in a ring-like shape. 
The CircularAnnulus function is used to define these ring-shaped areas. The annulus_aperture.to_mask function is used to mask these annular areas on the image data.

The sigma_clipped_stats function calculates the median value of the pixels in the annular area. 
This value is then used to calculate the aper_bkg (background value for the aperture image). 
The calc_total_error function is used to calculate the error value with the calculated background value.

The aperture_photometry function is used to calculate the light intensity inside the stars. This operation yields values such as aperture_sum and aperture_sum_err. 
These values are used to calculate the star's light intensity (magnitude) and error values, which are then written to a .csv file.

Finally, a matplotlib graph is created for each image, showing the star centers and annular areas. The "tb.close()" code is used to close the file. 
This process properly closes the file after the writing process, allowing for the release of open file resources when the operation is completed. 
This helps the program avoid unnecessary memory consumption and run more smoothly.
'''


tb=pd.read_csv('C:\\path\\tm.csv',sep=',')

JD=np.array(tb['JD'])-2459265.
fltr=((JD>0.4) & (JD<0.51)) 
JD=JD[fltr]
exptime=np.array(tb['exptime'])
x1=np.array(tb['x1'])
y1=np.array(tb['y1'])
x2=np.array(tb['x2'])
y2=np.array(tb['y2'])
mag1=np.array(tb['mag1'])
mag2=np.array(tb['mag2'])
emag1=np.array(tb['err1'])
emag2=np.array(tb['err2'])
dmag= mag2 - mag1

dmag=dmag[fltr] 
derr=np.sqrt(emag1**2.+emag2**2.) 
derr=derr[fltr] 
fig, ax = plt.subplots(1,1, figsize=(10,10)) 
ax.errorbar(JD, dmag,yerr=derr, ecolor='gray', elinewidth=0.5,marker='o',color='k',ls='none',ms=1)

poly_mod = PolynomialModel(1 ,prefix='p_') 
pars = poly_mod.guess(dmag, x=JD) 
gauss1 = GaussianModel(prefix='g1_') 
pars.update(gauss1.make_params()) 
pars['g1_center'].set(value=0.5, min=0.4, max=0.52) 
pars['g1_sigma'].set(value=0.01) 
pars['g1_height'].set(value=-0.8) 
mod = poly_mod + gauss1 

out = mod.fit(dmag, pars, x=JD)

print(out.fit_report(min_correl=0.5))

font1 = {'family': 'times new roman', 'color': 'black', 'size':'24'} 

plt.title('NY Vir',fontdict=font1)
plt.xlabel('Phase', fontdict=font1)
plt.ylabel('Normalized Diff. Mag.', fontdict=font1)
bp= pat.Patch(label='Gussian Fit')
ax.legend(handles=[bp])
ax.plot(JD,out.best_fit, lw=3,label='Gauss+Poly')
plt.show()


'''
----------------------------------TR----------------------------------
Bu kod, öncelikle bir veri seti olan 'tm.csv' dosyasını okuyarak verileri değişkenlere atar. Bu veri seti, bir yıldız çiftinin fotometrik verilerini içermektedir.

Verileri okuduktan sonra, verileri filtrelemek için bir dizi işlem yapılır. Öncelikle, tarihler JD'ye dönüştürülür ve verilerin yalnızca belirli bir aralığı seçilir. 
Bu filtreleme, belirli bir zaman aralığındaki verileri alarak ışık eğrisindeki belirli bir olayı analiz etmeyi amaçlar.

Daha sonra, değişkenleri belirleriz. x1, y1, x2, y2 değişkenleri, yıldızların pozisyonlarını verirken, mag1 ve mag2 değişkenleri, yıldızların parlaklıklarını ölçmek için kullanılan verileri içerir.

Sonrasında, iki yıldızın parlaklık değerleri arasındaki fark bulunur ve bu, ana yıldızın ışık eğrisi olarak kabul edilir. 
Daha sonra, bu veriler bir grafikte gösterilir ve polinom ve Gauss eğrisi modelleri oluşturmak için bir model örneği tanımlanır.

Model oluşturulduktan sonra, modellerle uyumlu parametreler tahmin edilir ve Gauss eğrisi için belirli sınırlar tanımlanır. Bu, modelleri gerçek verilere uygulamak için gereklidir.

Modellerin uygulanmasından sonra, grafik ve uyumlu verilerle birlikte en uygun model çizilir. Son olarak, grafik başlığı ve eksen etiketleri eklenir ve bir legend oluşturulur.

----------------------------------ENG----------------------------------
This code first reads a dataset, 'tm.csv', and assigns the data to variables. The dataset contains photometric data for a binary star system.

After reading the data, a series of operations are performed to filter the data. Firstly, the dates are converted to JD and only a certain range of data is selected. 
This filtering aims to analyze a specific event in the light curve by taking the data within a certain time range.

Next, variables are defined. The variables x1, y1, x2, and y2 give the positions of the stars, while mag1 and mag2 contain the data used to measure the brightness of the stars.

Then, the difference between the brightness values of the two stars is found, and this is considered as the light curve of the primary star. 
Subsequently, these data are plotted on a graph and a model instance is defined to create polynomial and Gaussian curve models. 

After the model is created, parameters that are consistent with the models are estimated, and specific boundaries are defined for the Gaussian curve. 
This is necessary to apply the models to the actual data.

After applying the models, the best model is plotted with the graph and the fit data. Finally, a title and axis labels are added to the graph, and a legend is created.
'''