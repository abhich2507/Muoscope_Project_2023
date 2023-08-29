#include <iostream>
#include <vector>
#include <cmath>

// Structure to represent a data point
struct DataPoint {
    double x;
    double y;
};

// Function to perform linear regression
void linearRegression(const std::vector<DataPoint>& data, double& slope, double& intercept) {
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXX = 0.0;
    int m= data.size();

    for (const DataPoint& point : data) {
        sumX += point.x;
        sumY += point.y;
        sumXY += point.x * point.y;
        sumXX += point.x * point.x;
    }

    double meanX = sumX / m;
    double meanY = sumY / m;

    slope = (sumXY - m* meanX * meanY) / (sumXX - m* meanX * meanX);
    intercept = meanY - slope * meanX;
}
