package ImageAnalysis;

import java.util.ArrayList;
import java.util.List;

import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfFloat;
import org.opencv.core.MatOfInt;
import org.opencv.core.Point;
import org.opencv.core.Scalar;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;

import ImageAnalysis.*;

public class HistogramCreator {

	List<Mat> istogrammi = new ArrayList<>();
	
	public List<Mat> getIstogrammi() { return istogrammi; }
	
	/** 
	 * crea gli istogrammi per ogni canale colore a partire dall'immagine LBP
	 * @param userWalk contiene le informazioni sulla camminata
	 * @param source e' il path da cui andare a caricare l'immagine LBP e in cui salvare le nuove immagini degli istogrammi
	 */
	public void createHistogramUserWalk(UserWalk userWalk, String source) {

		// carica la dll di opencv
		System.load("C:\\Users\\ggiorgi\\Desktop\\Tirocinio\\libs\\opencv\\build\\java\\x64\\opencv_java411.dll");
		
		// legge l'immagine LBP e la carica in un oggetto Mat di OpenCV
		Mat image = Imgcodecs.imread(source + userWalk.getWalk() + "lbp.png");

		// il ColorSpace standard id OpenCV e' BGR mentre le immagini LBP sono salvate in RGB, quindi bisogna applicare una conversione
		Imgproc.cvtColor(image, image, Imgproc.COLOR_BGR2RGB);
		
		// l'immagine LBP viene scomposta nei diversi canali colore
		List<Mat> rgbPlanes = new ArrayList<>();
		Core.split(image, rgbPlanes);

		// viene impostato il numero di bins dell'istogramma
		int histSize = 256;
		// viene impostato il range di valori. Il limite superiore e' esclusivo
		MatOfFloat histRange = new MatOfFloat(0, 256);

		// vengono istanziati gli oggetti che conterranno gli istogrammi
		Mat hist_r = new Mat();
		Mat hist_g = new Mat();
		Mat hist_b = new Mat();

		// vengono calcolati gli istogrammi tramite la funzione dedicata di OpenCV
		// i parametri di input sono la lista contenente le immagini da cui estrarre gli istogrammi, il canale colore, la maschera applicata (qui non usate, quindi l'oggetto creato e' vuoto)
		// l'oggetto Mat in cui andare a salvare l'istogramma, il numero di bins, il range di valori da usare, eventuali operazioni di aggregazione degli istogrammi (qui non usate)
		Imgproc.calcHist(rgbPlanes, new MatOfInt(0), new Mat(), hist_r, new MatOfInt(histSize), histRange, false);
		Imgproc.calcHist(rgbPlanes, new MatOfInt(1), new Mat(), hist_g, new MatOfInt(histSize), histRange, false);
		Imgproc.calcHist(rgbPlanes, new MatOfInt(2), new Mat(), hist_b, new MatOfInt(histSize), histRange, false);

		// vengono definiti i parametri dell'immagine in cui salvare l'istogramma
		int histW = 256, histH = 200;
		
		// vengono create le matrici delle immagini
		// i parametri di input sono l'altezza, la lunghezza, lo spazio colore e il colore di base
		// lo spazio colore qui usato � un generico composto da 3 canali (3C), ognuno di 8 bit unsigned (8U)
		Mat histImageR = new Mat(histH, histW, CvType.CV_8UC3, new Scalar(0, 0, 0));
		Mat histImageG = new Mat(histH, histW, CvType.CV_8UC3, new Scalar(0, 0, 0));
		Mat histImageB = new Mat(histH, histW, CvType.CV_8UC3, new Scalar(0, 0, 0));
		
		// gli istogrammi vengono normalizzati
		Core.normalize(hist_r, hist_r, 0, histImageR.rows(), Core.NORM_MINMAX);
		Core.normalize(hist_g, hist_g, 0, histImageG.rows(), Core.NORM_MINMAX);
		Core.normalize(hist_b, hist_b, 0, histImageB.rows(), Core.NORM_MINMAX);
		
		// i valori all'interno della matrice di ogni canale colore vengono riportati linearmente in un array dedicato 
		float[] rHistData = new float[(int) (hist_r.total() * hist_r.channels())];
		hist_r.get(0, 0, rHistData);
		float[] gHistData = new float[(int) (hist_g.total() * hist_g.channels())];
		hist_g.get(0, 0, gHistData);
		float[] bHistData = new float[(int) (hist_b.total() * hist_b.channels())];
		hist_b.get(0, 0, bHistData);
		
		// vengono disegnati gli istogrammi
		// i valori del colore per R e B sono invertiti, in quanto anche Scalar usa lo standard BGR
		for (int i = 0; i < histSize; i++) {
			//Imgproc.line(histImageB, new Point(i, histH), new Point(i, histH - Math.round(bHistData[i])), new Scalar(255, 0, 0), 1);
			
		}
		
		for (int i = 0; i < histSize; i++) {
		//	Imgproc.line(histImageG, new Point(i, histH), new Point(i, histH - Math.round(gHistData[i])), new Scalar(0, 255, 0), 1);
		//	Imgproc.line(histImageG, new Point(i, histH), new Point(i, histH - Math.round(gHistData[i])), new Scalar(0, 255, 0), 1);
		}
		
		for (int i = 0; i < histSize; i++) {
		//	Imgproc.line(histImageR, new Point(i, histH), new Point(i, histH - Math.round(rHistData[i])), new Scalar(0, 0, 255), 1);
		}
		
		
		// gli istogrammi vengono associati alla camminata. In questo modo vengono mantenuti in memoria e si velocizzano le operazioni di confronto
		List<Mat> istogrammi = new ArrayList<>();
		istogrammi.add(hist_r);
		istogrammi.add(hist_g);
		istogrammi.add(hist_b);
		
		userWalk.setHistograms(istogrammi);
		
		// vengono salvate le immagini degli istogrammi
		Imgcodecs.imwrite(source + "\\hist_r.png", histImageR);
		Imgcodecs.imwrite(source + "\\hist_g.png", histImageG);
		Imgcodecs.imwrite(source + "\\hist_b.png", histImageB);
	}
}
