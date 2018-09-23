package coreservlets;

import java.io.Serializable;

import javax.faces.bean.ManagedBean;

import org.primefaces.model.chart.CartesianChartModel;
import org.primefaces.model.chart.ChartSeries;

//@Named
@ManagedBean
public class ChartController implements Serializable {

    private CartesianChartModel model;

    public ChartController() {
        createLinearModel();
    }

    private void createLinearModel() {
        model = new CartesianChartModel();
        model.addSeries(getStockChartData("Stock Chart"));
    }

    private ChartSeries getStockChartData(String label) {
        ChartSeries data = new ChartSeries();
        data.setLabel(label);
        for (int i = 1; i <= 20; i++) {
            data.getData().put(i, (int) (Math.random() * 1000));
        }
        return data;
    }

    public CartesianChartModel getLinearModel() {
        return model;
    }

}