/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    if(theSource.getColumns() == 0 || theSource.getRows() == 0 || theTiles.empty())
    {
        return NULL;
    }

    MosaicCanvas * result = new MosaicCanvas(theSource.getRows(), theSource.getColumns());

    //I probably need to construct a KDTree with all my tiles
    vector<Point<3>> TileColors;
    map <Point<3>, TileImage*> TileMap; //in this way, we can map our point to our tile
    for(vector<TileImage>::iterator it = theTiles.begin(); it != theTiles.end(); ++it)
    {
        Point<3> averageColor = convertToXYZ(it->getAverageColor());
        TileMap[averageColor] = &(*it);
        TileColors.push_back(averageColor);
        
    }
    KDTree<3> colorTree(TileColors); //local variable, so I do not need to free memory on the heap

    //traverse the whole sourceImage with getRegion color and replace with nearestneighbor
    for(int i = 0; i < result->getRows(); i++)
    {
        for(int j = 0; j < result->getColumns(); j++)
        {
            Point<3> regionColor = convertToXYZ(theSource.getRegionColor(i,j));
            Point<3> neighbor = colorTree.findNearestNeighbor(regionColor);
            result->setTile(i, j, TileMap[neighbor]);
        }
    }


    return result;
}

