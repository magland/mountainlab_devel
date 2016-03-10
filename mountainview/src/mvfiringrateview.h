#ifndef MVFIRINGRATEVIEW_H
#define MVFIRINGRATEVIEW_H

#include <QWidget>
#include "diskreadmda.h"
#include "mvutils.h"

class MVFiringRateViewPrivate;
class MVFiringRateView : public QWidget
{
    Q_OBJECT
public:
    friend class MVFiringRateViewPrivate;
    MVFiringRateView();
    virtual ~MVFiringRateView();
    void setFirings(const DiskReadMda &X);
    void setVisibleLabels(const QList<int> &X);
    void setWindowWidth(double sec); //in seconds
    void setSamplingFreq(double ff); //in Hz
    void setCurrentEvent(MVEvent evt);
signals:

protected:
    //void paintEvent(QPaintEvent *evt);

private slots:
    void slot_update();

private:
    MVFiringRateViewPrivate *d;
};

#endif // MVFIRINGRATEVIEW_H
