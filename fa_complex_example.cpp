#include <bits/stdc++.h>
using namespace std;

#define XI 0.000001 // ξ

typedef double dbl;
typedef long long ll;
typedef pair<dbl, dbl> dd;
typedef vector<ll> vll;
typedef vector<dbl> vd;
typedef vector<dd> vdd;
typedef vector<vd> vvd;

const ll T = 100;			// number of iterations
const ll n = 20;			// number of locations (fireworks)
const ll d = 2;				// dimensionality
const ll m = 50;			// number of sparks with normal explosion
const ll mh = 50;			// number of sparks with Gaussian explosion
const dbl Ah = 20;			// maximum explosion amplitude
const dbl a = 0;		
const dbl b = 0.9;
const vdd bounds = { dd(-10, 10), dd(-10, 10) };	// Limits for d dimensions

// function d dimensions
dbl f(vd &v) { return - cos(v[0]) * cos(v[1]) * exp(- pow(v[0] - M_PI, 2) - pow(v[1] - M_PI, 2)); }

dbl randValue(dbl a, dbl b) 
{
	return ((dbl) rand() / (dbl) (RAND_MAX / (b - a))) + a;
}

dbl gaussianDistribution(dbl mean, dbl dev) 
{
	default_random_engine generator;
    generator.seed(time(NULL));
    normal_distribution<double> distribution(mean, dev);
    return distribution(generator);
}

dbl getDistance(vd &a, vd &b)
{
	dbl res = 0;
	for (int i = 0; i < a.size(); i++)
		res += pow(a[i] - b[i], 2);
	return sqrt(res);
}

void printSingle(vd &v) 
{
	cout << "[";
	for (int i = 0; i < v.size(); i++) 
	{
		cout << v[i];
		if (i < v.size() - 1) 
			cout << ", ";
	}
	cout << "]";
}

vll selectionDim(ll z)
{
	vll sel(d);
	for (int i = 0; i < sel.size(); i++)
	{
		if (randValue(0, 1) < 1.0 / (d - i)) {
			sel[i] = 1;
			if (!(--z)) break;
		}
	}
	return sel;
}

vvd generateLocations() 
{
	vvd locations;
	for (int i = 0; i < n; i++) 
	{
		vd loc;
		for (int j = 0; j < bounds.size(); j++) 
			loc.push_back(randValue(bounds[j].first, bounds[j].second));
		locations.push_back(loc);
	}
	return locations;
}

vd getFitness(vvd &fws, dbl &ymax, dbl &ymin, dbl &sYmaxF, dbl &sFYmin) 
{
	cout << "fireworks:" << endl;
	vd fitness;
	for (int i = 0; i < fws.size(); i++) 
	{
		dbl yfit = f(fws[i]);
		if (i == 0 or ymax < yfit)
			ymax = yfit;
		if (i == 0 or ymin > yfit)
			ymin = yfit;
		fitness.push_back(yfit);
		cout << "x_" << i+1 << " = "; printSingle(fws[i]); cout << " " << ", f(x_" << i+1 << ") = " << yfit << endl;
	}
	cout << "ymin = " << ymin << ", ymax = " << ymax << endl;
	sYmaxF = 0, sFYmin = 0;
	for (int i = 0; i < fitness.size(); i++)
	{
		sYmaxF += (ymax - fitness[i]);
		sFYmin += (fitness[i] - ymin);
	}
	return fitness;
}

ll numberSparksToFirework(dbl fit, dbl ymax, dbl sYmaxF)
{
	dbl s = m * ((ymax - fit + XI) / (sYmaxF + XI));
	s = (s < a * m ? a * m : (s > b * m ? b * m : s));
	return (ll) round(s);
}

dbl amplitudeOfExplotion(dbl fit, dbl ymin, dbl sFYmin)
{
	return Ah * ((fit - ymin + XI) / (sFYmin + XI));
}

vd getLocationSpark(vd &x, dbl A)
{
	vd xj = x;
	ll z = (ll) (round(d * randValue(0, 1)));
	vll sel = selectionDim(z);
	dbl h = A * randValue(-1, 1);
	for (int i = 0; i < xj.size(); i++)
	{
		if (sel[i]) 
		{
			xj[i] = xj[i] + h;
			if (xj[i] < bounds[i].first or xj[i] > bounds[i].second)
				xj[i] = bounds[i].first + (abs((ll)xj[i]) % ((ll)(bounds[i].second - bounds[i].first)));
		}
	}
	return xj;
}

void generateNormalSparks(vvd &spks, vd &x, ll sh, dbl A)
{
	for (int i = 0; i < sh; i++) 
	{
		vd spk = getLocationSpark(x, A);
		spks.push_back(spk);
		cout << "ns_"<< i+1 << " = "; printSingle(spk); cout << endl;
	}
}

void generateGaussianSparks(vvd &spks, vd &x)
{
	vd xj = x;
	ll z = (ll) (round(d * randValue(0, 1)));
	vll sel = selectionDim(z);
	dbl g = gaussianDistribution(1, 1);
	for (int i = 0; i < xj.size(); i++)
	{
		if (sel[i]) 
		{
			xj[i] = xj[i] * g;
			if (xj[i] < bounds[i].first or xj[i] > bounds[i].second)
				xj[i] = bounds[i].first + (abs((ll)xj[i]) % ((ll)(bounds[i].second - bounds[i].first)));
		}
	}
	spks.push_back(xj);
}

vd getBestLocation(vvd &locs)
{
	int best;
	for (int i = 0; i < locs.size(); i++)
	{
		if (i == 0 or f(locs[best]) > f(locs[i]))
			best = i;
	}
	return locs[best];
}

vvd calculateProbability(vvd &fws, vvd &nspks, vvd &gspks, vd &p)
{
	vvd joinAll = fws;
	joinAll.insert(joinAll.end(), nspks.begin(), nspks.end());
	joinAll.insert(joinAll.end(), gspks.begin(), gspks.end());
	
	vvd distances(joinAll.size(), vd(joinAll.size()));
	for (int i = 0; i < joinAll.size(); i++) {
		for (int j = i + 1; j < joinAll.size(); j++) {
			dbl dis = getDistance(joinAll[i], joinAll[j]);
			if (dis < 0)
				cout << "<--" << endl;
			distances[i][j] = dis;
			distances[j][i] = dis;
		}
	}

	dbl sumAllDistances = 0;
	vd sumDistances(joinAll.size());
	for (int i = 0; i < joinAll.size(); i++)
	{
		for (int j = 0; j < joinAll[i].size(); j++)
			sumDistances[i] += distances[i][j];
		sumAllDistances += sumDistances[i];
	}

	cout << endl << "probabilities:" << endl;
	p = vd(joinAll.size());
	for (int i = 0; i < joinAll.size(); i++)
	{
		dbl prob = sumDistances[i] / sumAllDistances;
		p[i] = (i == 0 ? prob : p[i - 1] + prob);
		cout << "l_" << i+1 << " = "; printSingle(joinAll[i]);
		cout << ", prob = " << prob << ", accum_prob = " << p[i] << endl;
	}
	return joinAll;
}

int main()
{
	srand(time(NULL));
	ios::sync_with_stdio(0); cout.tie(0);

	vd bestLocation;
	vvd locations = generateLocations();

	for (int i = 1; i <= T; i++) 
	{
		cout << "Iteration #" << i << endl << endl;
		dbl ymax, ymin, sumYmaxFunc, sumFuncYmin;
		vvd normalSparks, gaussianSparks;
		vvd fireworks = locations;
		vd fitness = getFitness(fireworks, ymax, ymin, sumYmaxFunc, sumFuncYmin);

		cout << endl << "normal sparks:" << endl;
		for (int j = 0; j < fireworks.size(); j++) 
		{
			ll sh = numberSparksToFirework(fitness[j], ymax, sumYmaxFunc);
			dbl A = amplitudeOfExplotion(fitness[j], ymin, sumFuncYmin);
			cout << "fws #" << j+1 << ": s = " << sh << ", A = " << A << endl;
			generateNormalSparks(normalSparks, fireworks[j], sh, A);
		}

		cout << endl << "gaussian sparks:" << endl;
		for (int j = 0; j < mh; j++) 
		{
			int r = rand() % n;
			generateGaussianSparks(gaussianSparks, fireworks[r]);
			cout << "gs_" << j+1 << " = "; printSingle(gaussianSparks[j]);
			cout << ", fw_" << r+1 << " = "; printSingle(fireworks[r]); cout << endl;
		}

		bestLocation = getBestLocation(locations);
		cout << endl << "Best Location:" << endl;
		cout << "x* = "; printSingle(bestLocation); cout << ", f(x*) = " << f(bestLocation) << endl;

		vd probabilitiesUnion;
		vvd unionFwsSpks = calculateProbability(fireworks, normalSparks, gaussianSparks, probabilitiesUnion);

		vvd newLocations(1, bestLocation);
		for (int j = 0; j < n - 1; j++)
		{
			int pos = lower_bound(probabilitiesUnion.begin(), probabilitiesUnion.end(), randValue(0, 1)) - probabilitiesUnion.begin();
			newLocations.push_back(unionFwsSpks[pos]);
		}
		locations = newLocations;

		cout << endl;
	}

	bestLocation = getBestLocation(locations);
	cout << "Best Location:" << endl;
	cout << "x* = "; printSingle(bestLocation); cout << ", f(x*) = " << f(bestLocation) << endl;

	return 0;
}

// Compile
// g++ -lm -O2 -std=c++11 fa_complex_example.cpp

// Execution
// ./a.out > complex_out.txt