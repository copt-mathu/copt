#include "PoissonDisk.h"
#include "gui.h"
#include <QApplication>

int main(int argc, char*argv[])
{
	// PoissonDisk pd(2,0.2,10);
	// Vector v1(2),v2(2);
	// v1[0]=-1;v1[1]=-1;
	// v2[0]=1;v2[1]=1;
	// pd.setRange(v1,v2);
	// pd.generate();
	QApplication a(argc,argv);
	QWidget *w = new PoissonWidget(30,20);
	w->show();
	a.exec();
}