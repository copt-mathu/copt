#include "gui.h"
#include "PoissonDisk.h"
#include <QPainter>
#include <QBrush>
#include <QGraphicsEllipseItem>


PoissonWidget::PoissonWidget(double r,int k,QWidget *parent)
	:
	QWidget(parent),
	__pd(2,r,k)
{
	__animation = new QGraphicsItemAnimation;
	__timeline = new TimeLine(1000,this);
	__timeline->setFrameRange(0,100);

	generatePoints();
}

QSize PoissonWidget::minimumSizeHint() const
{
	return QSize(500,500);
}

QSize PoissonWidget::sizeHint() const
{
	return QSize(500,500);
}

void PoissonWidget::keyPressEvent(QKeyEvent *e)
{
	if(e->key()==Qt::Key_W)
		generatePoints();
}

void PoissonWidget::generatePoints()
{
	Vector v1(2),v2(2);
	v1[0]=5;v1[1]=5;
	v2[0]=this->width()-5;v2[1]=this->height()-5;
	std::cout<<v1<<" "<<v2<<std::endl;
	__pd.setRange(v1,v2);
	__pd.generate();


	// update();
}

void PoissonWidget::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	painter.setRenderHint(QPainter::Antialiasing, true);
	painter.setPen(QPen(Qt::red));
	painter.setBrush(QBrush(Qt::red));
	for (auto p:__pd.points())
	{
		painter.drawEllipse(QPointF(p.coord[0],p.coord[1]),3,3);
		// std::cout<<p.coord<<std::endl;
	}
}