
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"
#include "stats.h"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
    :center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
    {}

toqutree::~toqutree(){
    clear(root);
}

toqutree::toqutree(const toqutree & other) {
    root = copy(other.root);
}


toqutree & toqutree::operator=(const toqutree & rhs){
    if (this != &rhs) {
        clear(root);
        root = copy(rhs.root);
    }
    return *this;
}

toqutree::toqutree(PNG & imIn, int k){

/* This constructor grabs the 2^k x 2^k sub-image centered */
/* in imIn and uses it to build a quadtree. It may assume  */
/* that imIn is large enough to contain an image of that size. */

    int sl = pow(2,k);
    PNG* centerIm = new PNG(sl,sl);
    int ulX = imIn.width()/2 - pow(2,k-1);
    int ulY = imIn.height()/2 - pow(2,k-1);

    for(int i = 0; i < sl; i++){
        for (int j = 0; j < sl; j++){
            HSLAPixel* p = imIn.getPixel(ulX+i, ulY+j);
            HSLAPixel* c = centerIm->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }
    root = buildTree(centerIm, k);
    delete(centerIm);
}

int toqutree::size() {
    return size(root);
}

int toqutree::size(const Node* croot){
    if(!croot){
        return 1;
    }
    return size(croot->NW) + size(croot->NE) + size(croot->SE) + size(croot->SW);
}


toqutree::Node * toqutree::buildTree(PNG * im, int k) {

/* your code here */

// Note that you will want to practice careful memory use
// In this function. We pass the dynamically allocated image
// via pointer so that it may be released after it is used .
// similarly, at each level of the tree you will want to
// declare a dynamically allocated stats object, and free it
// once you've used it to choose a split point, and calculate
// an average.

    stats* im_stats = new stats(*im);
    Node* croot;
    int sl = pow(2,k);
    if(k == 0){
        croot = new Node(make_pair(0,0), 0, *im->getPixel(0,0));
        croot->SE = NULL;
        croot->NE = NULL;
        croot->NW = NULL;
        croot->SW = NULL;
        delete(im_stats);
        return croot;
    }

    if(k == 1) {
        pair<int, int> ul = make_pair(0, 0);
        pair<int, int> lr = make_pair(1, 1);

        croot = new Node(lr, 1, im_stats->getAvg(ul, lr));

        PNG* SE = new PNG(1,1);
		PNG* SW = new PNG(1,1);
		PNG* NE = new PNG(1,1);
		PNG* NW = new PNG(1,1);

		HSLAPixel* sePixel = im->getPixel(1,1);
		HSLAPixel* swPixel = im->getPixel(0,1);
		HSLAPixel* nePixel = im->getPixel(1,0);
		HSLAPixel* nwPixel = im->getPixel(0,0);

		HSLAPixel* se = SE->getPixel(0,0);
		HSLAPixel* sw = SW->getPixel(0,0);
		HSLAPixel* ne = NE->getPixel(0,0);
		HSLAPixel* nw = NW->getPixel(0,0);

		*se = *sePixel;
		*sw = *swPixel;
		*ne = *nePixel;
		*nw = *nwPixel;

		croot->SE = buildTree(SE, 0);
		croot->SW = buildTree(SW, 0);
		croot->NE = buildTree(NE, 0);
		croot->NW = buildTree(NW, 0);
        delete(im_stats);
        delete(SE);
        delete(SW);
        delete(NE);
        delete(NW);
        return croot;
    }

    double min_entropy = 100000000000;
    pair<int, int> split_point = make_pair(0, 0);

    // case 1 - traverse through points in Quadrant II
    for(int i = sl/4; i < sl/2; i++){
        for (int j = sl/4; j < sl/2; j++){
            vector<int> SE_Hist = im_stats->buildHist(make_pair(i,j),make_pair(i+sl/2-1,j+sl/2-1));
            vector<int> NE_Hist0 = im_stats->buildHist(make_pair(i,j+sl/2), make_pair(i+sl/2-1,sl-1));
            vector<int> NE_Hist1 = im_stats->buildHist(make_pair(i,0),make_pair(i+sl/2-1,j-1));
            vector<int> SW_Hist0 = im_stats->buildHist(make_pair(i+sl/2,j),make_pair(sl-1,j+sl/2-1));
            vector<int> SW_Hist1 = im_stats->buildHist(make_pair(0,j),make_pair(i-1, j+sl/2-1));
            vector<int> NW_Hist0 = im_stats->buildHist(make_pair(i+sl/2,j+sl/2),make_pair(sl-1,sl-1));
            vector<int> NW_Hist1 = im_stats->buildHist(make_pair(0,j+sl/2),make_pair(i-1,sl-1));
            vector<int> NW_Hist2 = im_stats->buildHist(make_pair(0,0), make_pair(i-1,j-1));
            vector<int> NW_Hist3 = im_stats->buildHist(make_pair(i+sl/2,0), make_pair(sl-1, j-1));
            vector<int> NE_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = NE_Hist0[i] + NE_Hist1[i];
                NE_Hist.push_back(sum);
            }
            vector<int> SW_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = SW_Hist0[i] + SW_Hist1[i];
                SW_Hist.push_back(sum);
            }
            vector<int> NW_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = NW_Hist0[i] + NW_Hist1[i] + NW_Hist2[i] + NW_Hist3[i];
                NW_Hist.push_back(sum);
            }
            int area = (pow(2,k-1) * pow(2,k-1));
            double SE_entropy = im_stats->entropy(SE_Hist, area);
            double NE_entropy = im_stats->entropy(NE_Hist, area);
            double SW_entropy = im_stats->entropy(SW_Hist, area);
            double NW_entropy = im_stats->entropy(NW_Hist, area);

            double total_entropy = SE_entropy + NE_entropy + SW_entropy + NW_entropy;
            if (total_entropy < min_entropy){
                min_entropy = total_entropy;
                split_point = make_pair(i, j);
            }
        }
    }

    //case 2 - traverse through points in Quadrant I
    int Q1Shift = 0;
    for(int i = sl/2; i < 3*sl/4; i++){
        Q1Shift = i - sl/2;
        for (int j = sl/4; j < sl/2; j++){
            vector<int> SE_Hist0 = im_stats->buildHist(make_pair(i,j), make_pair(sl-1,j+sl/2-1));
            vector<int> SE_Hist1 = im_stats->buildHist(make_pair(0,j), make_pair(Q1Shift-1, j+sl/2-1));
            vector<int> SW_Hist = im_stats->buildHist(make_pair(i-sl/2,j),make_pair(i-1,j+sl/2 -1));
            vector<int> NE_Hist0 = im_stats->buildHist(make_pair(i,0), make_pair(sl-1,j-1));
            vector<int> NE_Hist1 = im_stats->buildHist(make_pair(i,j+sl/2),make_pair(sl-1, sl-1));
            vector<int> NE_Hist2 = im_stats->buildHist(make_pair(0,0),make_pair(Q1Shift-1,j-1));
            vector<int> NE_Hist3 = im_stats->buildHist(make_pair(0,j+sl/2),make_pair(Q1Shift-1, sl-1));
            vector<int> NW_Hist0 = im_stats->buildHist(make_pair(i-sl/2,0),make_pair(i-1,j-1));
            vector<int> NW_Hist1 = im_stats->buildHist(make_pair(i-sl/2,j+sl/2),make_pair(i-1,sl-1));
            vector<int> SE_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = SE_Hist0[i] + SE_Hist1[i];
                SE_Hist.push_back(sum);
            }
            vector<int> NW_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = NW_Hist0[i] + NW_Hist1[i];
                NW_Hist.push_back(sum);
            }
            vector<int> NE_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = NE_Hist0[i] + NE_Hist1[i] + NE_Hist2[i] + NE_Hist3[i];
                NE_Hist.push_back(sum);
            }
            int area = (pow(2,k-1) * pow(2,k-1));
            double SE_entropy = im_stats->entropy(SE_Hist, area);
            double NE_entropy = im_stats->entropy(NE_Hist, area);
            double SW_entropy = im_stats->entropy(SW_Hist, area);
            double NW_entropy = im_stats->entropy(NW_Hist, area);
            double total_entropy = SE_entropy + NE_entropy + SW_entropy + NW_entropy;
            if (total_entropy < min_entropy){
                min_entropy = total_entropy;
                split_point = make_pair(i, j);
            }
        }
    }

    // case 3 - traverse throught points in Quadrant III
    int Q3Shift = 0;
    for (int j = sl/2; j < 3*sl/4; j++){
        for(int i = sl/4; i < sl/2; i++){
            vector<int> SE_Hist0 = im_stats->buildHist(make_pair(i,j),make_pair(i+sl/2 -1, sl-1));
            vector<int> SE_Hist1 = im_stats->buildHist(make_pair(i,0),make_pair(i+sl/2 - 1, Q3Shift-1));
            vector<int> SW_Hist0 = im_stats->buildHist(make_pair(0,0),make_pair(i-1,Q3Shift-1));
            vector<int> SW_Hist1 = im_stats->buildHist(make_pair(i+sl/2,0),make_pair(sl-1,Q3Shift-1));
            vector<int> SW_Hist2 = im_stats->buildHist(make_pair(0,j),make_pair(i-1,sl-1));
            vector<int> SW_Hist3 = im_stats->buildHist(make_pair(i+sl/2,j),make_pair(sl-1,sl-1));
            vector<int> NE_Hist = im_stats->buildHist(make_pair(i, j-sl/2),make_pair(i+sl/2-1, j-1));
            vector<int> NW_Hist0 = im_stats->buildHist(make_pair(i+sl/2, j-sl/2),make_pair(sl-1,j-1));
            vector<int> NW_Hist1 = im_stats->buildHist(make_pair(0, j-sl/2),make_pair(i-1,j-1));
            vector<int> SE_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = SE_Hist0[i] + SE_Hist1[i];
                SE_Hist.push_back(sum);
            }
            vector<int> SW_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = SW_Hist0[i] + SW_Hist1[i] + SW_Hist2[i] + SW_Hist3[i];
                SW_Hist.push_back(sum);
            }
            vector<int> NW_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = NW_Hist0[i] + NW_Hist1[i];
                NW_Hist.push_back(sum);
            }
            int area = (pow(2,k-1) * pow(2,k-1));
            double SE_entropy = im_stats->entropy(SE_Hist, area);
            double NE_entropy = im_stats->entropy(NE_Hist, area);
            double SW_entropy = im_stats->entropy(SW_Hist, area);
            double NW_entropy = im_stats->entropy(NW_Hist, area);
            double total_entropy = SE_entropy + NE_entropy + SW_entropy + NW_entropy;
            if (total_entropy < min_entropy){
                min_entropy = total_entropy;
                split_point = make_pair(i, j);
            }
        }
        Q3Shift++;
    }


    //mcase 4 - traverse throught points in Quadrant IV
    int Q4XShift = 0;
    for (int i = sl/2; i < 3*sl/4; i++){
        int Q4YShift = 0; 
        for(int j = sl/2; j < 3*sl/4; j++){
            vector<int> SE_Hist0 = im_stats->buildHist(make_pair(i,j),make_pair(sl-1,sl-1));
            vector<int> SE_Hist1 = im_stats->buildHist(make_pair(i,0),make_pair(sl-1, Q4YShift-1));
            vector<int> SE_Hist2 = im_stats->buildHist(make_pair(0,0),make_pair(Q4XShift-1,Q4YShift-1));
            vector<int> SE_Hist3 = im_stats->buildHist(make_pair(0,j),make_pair(Q4XShift-1, sl-1));
            vector<int> SW_Hist0 = im_stats->buildHist(make_pair(i-sl/2,j),make_pair(i-1,sl-1));
            vector<int> SW_Hist1 = im_stats->buildHist(make_pair(i-sl/2,0),make_pair(i-1, Q4YShift-1));
            vector<int> NE_Hist0 = im_stats->buildHist(make_pair(i,j-sl/2),make_pair(sl-1,j-1));
            vector<int> NE_Hist1 = im_stats->buildHist(make_pair(0, j-sl/2),make_pair(Q4XShift-1,j-1));
            vector<int> NW_Hist = im_stats->buildHist(make_pair(i-sl/2,j-sl/2),make_pair(i-1,j-1));
            vector<int> NE_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = NE_Hist0[i] + NE_Hist1[i];
                NE_Hist.push_back(sum);
            }
            vector<int> SW_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = SW_Hist0[i] + SW_Hist1[i];
                SW_Hist.push_back(sum);
            }
            vector<int> SE_Hist;
            for (int i = 0; i < 36 ; i++) {
                int sum = SE_Hist0[i] + SE_Hist1[i] + SE_Hist2[i] + SE_Hist3[i];
                SE_Hist.push_back(sum);
            }
            int area = (pow(2,k-1) * pow(2,k-1));
            double SE_entropy = im_stats->entropy(SE_Hist, area);
            double NE_entropy = im_stats->entropy(NE_Hist, area);
            double SW_entropy = im_stats->entropy(SW_Hist, area);
            double NW_entropy = im_stats->entropy(NW_Hist, area);
            double total_entropy = SE_entropy + NE_entropy + SW_entropy + NW_entropy;
            if (total_entropy < min_entropy){
                min_entropy = total_entropy;
                split_point = make_pair(i, j);
            }
            Q4YShift++;
        }
        Q4XShift++;
    }

    int newL = pow(2,k-1);
    int SE_X = split_point.first;
    int SE_Y = split_point.second;

    // set up SE 
    PNG* SE_PNG = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+i)%sl;
            int y = (SE_Y+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = SE_PNG->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }
    Node* SE_child = buildTree(SE_PNG, k-1);

    // set up NE
    PNG* NE_PNG = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+i)%sl;
            int y = (SE_Y+newL+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = NE_PNG->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }
    Node* NE_child = buildTree(NE_PNG, k-1);
    

    // set up SW
    PNG* SW_PNG = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+newL+i)%sl;
            int y = (SE_Y+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = SW_PNG->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }
    Node* SW_child = buildTree(SW_PNG, k-1);

    // set up NW
    PNG* NW_PNG = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+newL+i)%sl;
            int y = (SE_Y+newL+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = NW_PNG->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }
    Node* NW_child = buildTree(NW_PNG, k-1);

    croot = new Node(split_point,k,im_stats->getAvg(make_pair(0,0), make_pair((int)im->width() - 1, (int)im->height() - 1)));

    croot->SE = SE_child;
    croot->NW = NW_child;
    croot->NE = NE_child;
    croot->SW = SW_child;
    delete(im_stats);
    delete(SE_PNG);
    delete(SW_PNG);
    delete(NE_PNG);
    delete(NW_PNG);
    return croot;
}


PNG toqutree::render(){

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the
// quadtree, instead.

/* your code here */
   	PNG png (pow(2, root->dimension), pow(2, root->dimension));
	for (int i = 0; i < (int)png.width(); i++){
		for (int j = 0; j < (int)png.height(); j++){
			pair<int, int> cord(i, j);
			Node* node = find(root, cord);
			*png.getPixel(i, j) = node->avg;
		}
	}
	return png;
}

toqutree::Node * toqutree::find(Node* croot, const pair<int, int> &key){
    // todo
    // croot is leaf
    pair<int, int> ctr = croot->center;
	int k = croot->dimension;
    int origRend = pow(2,k) - 1;
    int subLength = pow(2, k - 1);
	if (k == 0 || croot->SE == NULL) return croot;

    if(k == 1){
        if (key.first == 0){
            if(key.second == 0 ){
                return find(croot->NW, make_pair(0,0));
            }else{
                return find(croot->SW, make_pair(0,0));
            }
        }else{
            if(key.second == 0){
                return find(croot->NE, make_pair(0,0));
            }else{
                return find(croot->SE, make_pair(0,0));
            }
        }
    }

    // SE IN QI
	else if(((ctr.first + subLength) > origRend) && ((ctr.second + subLength) <=  origRend)){
		if((((key.first >= ctr.first) && (key.first <= pow(2, k) - 1)) || ((key.first >= 0) && (key.first < ctr.first-subLength))) && 
            (key.second >= ctr.second) && (key.second < ctr.second+subLength)){
                int newX = ctr.first;
                int newY = ctr.second;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(croot->SE, newKey);
            }
		else if((((key.second >= ctr.second + subLength) && (key.second <= pow(2, k) - 1)) || ((key.second>=0) && (key.second < ctr.second))) &&
            (key.first >= ctr.first - subLength) && (key.first < ctr.first)){
                int newX = ctr.first - subLength;
                int newY = ctr.second + subLength;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(croot->NW, newKey);
            }
		else if((key.second >= ctr.second) && (key.second < (ctr.second + subLength)) && (key.first < ctr.first) && (key.first >= (ctr.first - subLength))){
            int newX = ctr.first - subLength;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->SW, newKey);
        }
		else{
            int newX = ctr.first;
            int newY = ctr.second +subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->NE, newKey);
        }
			
	}

    // SE IN QII
	else if(((ctr.first + subLength) <= origRend) && ((ctr.second + subLength) <= origRend)){
		if((key.first < (ctr.first + subLength)) && (key.second < (ctr.second + subLength)) && (key.first >= ctr.first) && (key.second >= ctr.second)){
            pair<int,int> newKey = make_pair(key.first-ctr.first, key.second-ctr.second);
            return find(croot->SE, newKey);
        }
		else if((key.first >= ctr.first) && (key.first < (ctr.first + subLength)) && 
        (((key.second < ctr.second) && (key.second >=0)) || ((key.second< pow(2, k)) && (key.second >= (ctr.second + subLength))))){
			int newX = ctr.first;
            int newY = ctr.second+subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->NE, newKey);
        }
		else if((key.second >= ctr.second) && (key.second < (ctr.second + subLength)) && ((key.first < ctr.first) || (key.first >= (ctr.first + subLength)))){
			int newX = ctr.first + subLength;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->SW, newKey);
        }
		else{
            int newX = ctr.first +subLength;
            int newY = ctr.second + subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->NW, newKey);
        }
	}


    // SE IN QIII
	else if(((ctr.first + subLength) <= origRend) && ((ctr.second + subLength) > pow(2, k) - 1)){
        if((key.first < (ctr.first + subLength)) && (key.first >= ctr.first) && 
            ((key.second >= ctr.second && key.second <= pow(2, k) - 1) || (key.second>=0 && key.second < (ctr.second - subLength)))){
                int newX = ctr.first;
                int newY = ctr.second;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(croot->SE, newKey);
            }
        else if((key.first >= ctr.first) && (key.first < (ctr.first + subLength)) && (key.second < ctr.second) && (key.second >= (ctr.second - subLength))){
                int newX = ctr.first;
                int newY = ctr.second -subLength;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(croot->NE, newKey);
        }
        else if((key.second < ctr.second) && (key.second >= (ctr.second - subLength)) && 
            (((key.first >= ctr.first+subLength) && (key.first < pow(2, k))) || ((key.first>=0) && (key.first < ctr.first)))){
                int newX = ctr.first + subLength;
                int newY = ctr.second -subLength;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(croot->NW, newKey);
            }
		else{
                int newX = ctr.first + subLength;
                int newY = ctr.second;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(croot->SW, newKey);
        }
	}
    
    // SE IN QIV
    else {
		if((key.first >= (ctr.first - subLength)) && (key.second >= (ctr.second - subLength)) && (key.first < ctr.first) && (key.second < ctr.second)){
            int newX = ctr.first - subLength;
            int newY = ctr.second - subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->NW, newKey);
        }
		else if((((key.first >= ctr.first) && (key.first < pow(2, k))) || ((key.first >= 0) && (key.first < ctr.first-subLength))) && 
        (key.second >= ctr.second - subLength) && (key.second < ctr.second)){
            int newX = ctr.first;
            int newY = ctr.second -subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->NE, newKey);
        }
		else if((key.first >= (ctr.first - subLength)) && (key.first < ctr.first)){
            int newX = ctr.first-subLength;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->SW, newKey);
        }
		else{
            int newX = ctr.first;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(croot->SE, newKey);
        }
	}
}

/* oops, i left the implementation of this one in the file! */
void toqutree::prune(double tol){
    prune(root,tol);
}

void toqutree::prune(Node * croot, double tol){
    if(croot != NULL){
        HSLAPixel avg = croot->avg;
        if(canPrune(croot, tol, avg)){
            clear(croot->SE);
            clear(croot->SW);
            clear(croot->NE);
            clear(croot->NW);
        } else{
            prune(croot->SW, tol);
            prune(croot->NE, tol);
            prune(croot->NW, tol);
            prune(croot->SE, tol);
        }
        
    }

}

bool toqutree::canPrune(Node* croot, double tol, HSLAPixel avgPixel){
    if(croot->SE == NULL  && croot->SW == NULL && croot->NE == NULL && croot->NW == NULL ){
        if(avgPixel.dist(croot->avg) <= tol){
            return true;
        } else{
            return false;
        }
    }else{
            return canPrune(croot->SE, tol, avgPixel)
                && canPrune(croot->SW, tol, avgPixel) 
                && canPrune(croot->NE, tol, avgPixel)
                && canPrune(croot->NW, tol, avgPixel);
    }
}




/* called by destructor and assignment operator*/
void toqutree::clear(Node * & curr){
    if(curr == NULL) return;
    clear(curr->NW);
    clear(curr->NE);
    clear(curr->SE);
    clear(curr->SW);
    delete(curr);
    curr = NULL;
}

/* done */
/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {
    Node* root = NULL;
    if (other != NULL) {
        root = new Node(other->center, other->dimension, other->avg);
        root->SE = copy(other->SE);
		root->SW = copy(other->SW);
		root->NE = copy(other->NE);
		root->NW = copy(other->NW);
    }
    return root;
}


pair<int,int> toqutree::findNewKey(pair<int,int> key, int k, int X, int Y){
    int newF;
    int newS;
    
    if(key.first >= X){
        newF = key.first - X;
    }else{
        newF = key.first + (pow(2,k) - X);
    }

    if(key.second >= Y){
        newS = key.second - Y;
    }else{
        newS = key.second + (pow(2,k) - Y);
    }

    return make_pair(newF, newS);
}