#ifndef POISSON_GUI_H__
#define POISSON_GUI_H__

#include <QWidget>
#include "PoissonDisk.h"
#include <QKeyEvent>
#include <QTimeLine>
#include <QGraphicsItemAnimation>

class PoissonWidget
	:
	public QWidget
{
	Q_OBJECT
private:

	PoissonDisk 				__pd;
	QTimeLine					*__timeline;
	QGraphicsItemAnimation 		*__animation;
public:
	PoissonWidget(double r,int k,QWidget *parent = 0);
	~PoissonWidget(){}

	QSize minimumSizeHint() const;
	QSize sizeHint() const;

	
	void keyPressEvent(QKeyEvent *e);
	void paintEvent(QPaintEvent*);

signals:

public slots:
	void generatePoints();
};

#endif