
#include "stats.h"

stats::stats(PNG & im){
     int height = im.height();
     int width = im.width();
     // initialize size (resize vectors)
     for (int i = 0; i < width; i++){
         vector<double> columnS;
         vector<double> columnL;
         vector<double> columnHueX;
         vector<double> columnHueY;
         vector<vector<int>> columnHist;
         for (int j = 0; j <height; j++){
            columnS.push_back(0);
            columnL.push_back(0);
            columnHueX.push_back(0);
            columnHueY.push_back(0);
            vector<int> freq;
            for (int k = 0; k < 36; k++){
                freq.push_back(0);
            }
            columnHist.push_back(freq);
         }
         sumSat.push_back(columnS);
         sumLum.push_back(columnL);
         sumHueX.push_back(columnHueX);
         sumHueY.push_back(columnHueY);
         hist.push_back(columnHist);
    }
    // calculate position (0,0) 
     sumSat[0][0] = im.getPixel(0,0)->s;
     sumLum[0][0] = im.getPixel(0,0)->l;
     sumHueX[0][0] = cos(im.getPixel(0,0)->h * PI / 180);
     sumHueY[0][0] = sin(im.getPixel(0,0)->h * PI / 180);
     int h = (int) (im.getPixel(0,0)->h / 10);
     hist[0][0][h] = 1;
    // calculate first row
     for (int i = 1; i < width; i++){
         sumSat[i][0] = sumSat[i-1][0] + im.getPixel(i,0)->s;
         sumLum[i][0] = sumLum[i-1][0] + im.getPixel(i,0)->l;
         sumHueX[i][0] = sumHueX[i-1][0] + cos(im.getPixel(i,0)->h * PI / 180);
         sumHueY[i][0] = sumHueY[i-1][0] + sin(im.getPixel(i,0)->h * PI / 180);
         for (int k = 0; k < 36; k++){
             hist[i][0][k] = hist[i-1][0][k];
         }
         h = (int) (im.getPixel(i,0)->h /10);
         hist[i][0][h]++;  
     }
     // calculate first coloumn
     for (int j = 1; j < height; j++){
         sumSat[0][j] = sumSat[0][j-1] + im.getPixel(0,j)->s;
         sumLum[0][j] = sumLum[0][j-1] + im.getPixel(0,j)->l;
         sumHueX[0][j] = sumHueX[0][j-1] + cos(im.getPixel(0,j)->h * PI / 180);
         sumHueY[0][j] = sumHueY[0][j-1] + sin(im.getPixel(0,j)->h * PI / 180);
         for (int k = 0; k < 36; k++){
             hist[0][j][k] = hist[0][j-1][k];
         }
         h = (int) (im.getPixel(0,j)->h /10);
         hist[0][j][h]++;  
     }
     // calculate rest
     for (int i = 1; i < width; i++){
         for (int j = 1; j <height; j++){
             sumSat[i][j] = sumSat[i-1][j] + sumSat[i][j-1] - sumSat[i-1][j-1] + im.getPixel(i,j)->s;
             sumLum[i][j] = sumLum[i-1][j] + sumLum[i][j-1] - sumLum[i-1][j-1] + im.getPixel(i,j)->l;
             sumHueX[i][j] = sumHueX[i-1][j] + sumHueX[i][j-1] - sumHueX[i-1][j-1] + cos(im.getPixel(i,j)->h / 180 * PI);
			 sumHueY[i][j] = sumHueY[i-1][j] + sumHueY[i][j-1] - sumHueY[i-1][j-1] + sin(im.getPixel(i,j)->h / 180 * PI);
            for(int k = 0; k < 36; k++){
                 hist[i][j][k] = hist[i-1][j][k] + hist[i][j-1][k] - hist[i-1][j-1][k];
            }
            h = (int) (im.getPixel(i,j)->h / 10);
            hist[i][j][h]++;
         }
     }
}

long stats::rectArea(pair<int,int> ul, pair<int,int> lr){
    int upperX = ul.first;
    int upperY = ul.second;
    int lowerX = lr.first;
    int lowerY = lr.second;
    return (long) (lowerX - upperX + 1) * (lowerY - upperY + 1);
}

HSLAPixel stats::getAvg(pair<int,int> ul, pair<int,int> lr){
/* your code here */
    long area = rectArea(ul,lr);
    HSLAPixel avgPixel;
    double sumHX = sumHueX[lr.first][lr.second];
    double sumHY = sumHueY[lr.first][lr.second];
    double sumS = sumSat[lr.first][lr.second];
    double sumL = sumLum[lr.first][lr.second];
    if(ul.first){
        sumHX -= sumHueX[ul.first - 1][lr.second];
        sumHY -= sumHueY[ul.first - 1][lr.second];
        sumS -= sumSat[ul.first - 1][lr.second];
        sumL -= sumLum[ul.first - 1][lr.second];
    }
    if(ul.second){
        sumHX -= sumHueX[lr.first][ul.second - 1];
        sumHY -= sumHueY[lr.first][ul.second - 1];
        sumS -= sumSat[lr.first][ul.second - 1];
        sumL -= sumLum[lr.first][ul.second - 1];
        
        if(ul.first){
            sumHX += sumHueX[ul.first - 1][ul.second - 1];
            sumHY += sumHueY[ul.first - 1][ul.second - 1];
            sumS += sumSat[ul.first - 1][ul.second - 1];
            sumL += sumLum[ul.first - 1][ul.second - 1];
        }
    }
    double avgHX = (double) sumHX/area;
    double avgHY = (double) sumHY/area;
    double avgH = atan2(avgHY, avgHX) * 180 / PI;
    if(avgH <0){
        avgH += 360;
    }
    double avgS = (double) sumS/area;
    double avgL = (double) sumL/area;
    avgPixel.a = 1.0;
    avgPixel.s = avgS;
    avgPixel.l = avgL;
    avgPixel.h = avgH;
    return avgPixel;
}

vector<int> stats::buildHist(pair<int,int> ul, pair<int,int> lr){
    vector<int> distn;
    for(int i = 0; i < 36; i++){
            distn.push_back(0);
        }
    if(lr.second < ul.second || lr.first<ul.first || lr.second<0 || lr.first<0 || ul.first<0 || ul.second<0){
        return distn;
    }
    for (int i = 0; i < 36; i++){
        distn[i] = hist[lr.first][lr.second][i];
    }
    if (ul.first != 0){
        for (int i = 0; i < 36; i++){
            distn[i] = distn[i] - hist[ul.first-1][lr.second][i];
        }
    }
    if(ul.second != 0){
        for (int i = 0; i < 36; i++){
            distn[i] = distn[i] - hist[lr.first][ul.second-1][i];
        }
        if (ul.first) {
			for (int i = 0; i < 36; i++)
				distn[i] = distn[i] + hist[ul.first - 1][ul.second - 1][i];
		}
    }
    return distn;
/* your code here */
}

// takes a distribution and returns entropy
// partially implemented so as to avoid rounding issues.
double stats::entropy(vector<int> & distn,int area){

    double entropy = 0.;

/* your code here */

    for (int i = 0; i < 36; i++) {
        if (distn[i] > 0 ) 
            entropy += ((double) distn[i]/(double) area) 
                                    * log2((double) distn[i]/(double) area);
    }

    return  -1 * entropy;

}

double stats::entropy(pair<int,int> ul, pair<int,int> lr){
    vector<int> distn = buildHist(ul,lr);
    return entropy(distn,rectArea(ul, lr));
}
