import { Grid } from "@mantine/core";
import { Reaction } from "@/lib/types";
import {
  MonomerInfoCard,
  ProductInfoCard,
  ParametersInfoCard,
  ThermoInfoCard,
  MetadataInfoCard,
} from "@/components/InfoCards";
import {
  VanthoffPlotCard,
  ExtrapolationPlotCard,
} from "@/components/PlotCards";

type Props = {
  data: Reaction;
};

export default function ReactionSummaryPage({ data }: Props) {
  const hasVanthoff = data.thermo?.vanthoff;
  const hasExtrapolation = data.thermo?.extrapolation;

  return (
    <Grid gutter="xl">
      <Grid.Col span={6}>
        <MonomerInfoCard data={data.monomer} />
        <ParametersInfoCard data={data.parameters} />
      </Grid.Col>
      <Grid.Col span={6}>
        <ProductInfoCard data={data.product} />
        <ThermoInfoCard data={data.thermo} />
        <MetadataInfoCard data={data.metadata} />
      </Grid.Col>

      {hasVanthoff && (
        <Grid.Col span={12}>
          <VanthoffPlotCard data={data.thermo.vanthoff} />
        </Grid.Col>
      )}

      {hasExtrapolation && (
        <Grid.Col span={12}>
          <ExtrapolationPlotCard data={data.thermo.extrapolation} />
        </Grid.Col>
      )}
    </Grid>
  );
}
