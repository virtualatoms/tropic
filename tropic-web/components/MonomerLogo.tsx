import { Card, Badge, Stack, Box, Center, Image } from "@mantine/core";

export function MonomerLogo({
  svg,
  monomerId,
}: {
  svg: string;
  monomerId: string;
}) {
  return (
    <Box mb={20} px={0}>
      <Stack>
        <Card
          withBorder
          radius="lg"
          style={{
            borderColor: "var(--mantine-color-blue-5)",
            borderWidth: 0,
          }}
        >
          <Image src={`data:image/svg+xml;base64,${svg}`} alt={monomerId} />
        </Card>
        <Center style={{ marginTop: -10, zIndex: 1 }}>
          <Badge size="lg" radius="lg" variant="outline">
            {monomerId}
          </Badge>
        </Center>
      </Stack>
    </Box>
  );
}
