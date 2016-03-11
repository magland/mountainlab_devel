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
    void setTimesAmplitudes(const QList<double> &times,const QList<double> &amps);
    void setSamplingFreq(double ff); //in Hz
    void setCurrentEvent(MVEvent evt);
    void setEpochs(const QList<Epoch> &epochs);
signals:

protected:
    void paintEvent(QPaintEvent *evt);

private slots:
    void slot_update();

private:
    MVFiringRateViewPrivate *d;
};

#endif // MVFIRINGRATEVIEW_H
