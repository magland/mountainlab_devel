#ifndef MVCLIPSVIEW_H
#define MVCLIPSVIEW_H

#include "sstimeseriesview.h"

class MVClipsViewPrivate;
class MVClipsView : public SSTimeSeriesView
{
    Q_OBJECT
public:
    friend class MVClipsViewPrivate;
    MVClipsView();
    virtual ~MVClipsView();
    void setClips(const Mda &clips);
    void setTimes(const QList<double> &times);
    int currentClipIndex();
    double currentClipTimepoint();
signals:
    void currentClipTimepointChanged();
public slots:
private:
    MVClipsViewPrivate *d;
};

#endif // MVCLIPSVIEW_H
